import os
import tempfile
import shutil
import multiprocessing
from typing import Optional, List

import rdkit.Chem as Chem
from pydantic import BaseModel
from typing_extensions import Literal

from dockstream.core.docker import Docker
from dockstream.core.Dock37.Dock37_result_parser import Dock37ResultParser
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.execute_external.Dock37 import Dock37Executor
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.dockstream_exceptions import DockingRunFailed

_LP = RDkitLigandPreparationEnum()
_LE = LoggingConfigEnum()

class Dock37Parameters(BaseModel):
    check_dir: Optional[str] = None
    dockfiles_indock_dir: Optional[List[str]] = None
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None


class Dock37(Docker, BaseModel):
    """Interface to the "UCSF DOCK3.7" backend."""

    backend: Literal["Dock37"] = "Dock37"
    parameters: Dock37Parameters

    _DOCK37_executor: Dock37Executor = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

    def _initialize_executors(self):
        """Initialize executors and check if they are available."""

        self._DOCK37_executor = Dock37Executor(
            prefix_execution=self.parameters.prefix_execution, 
            binary_location=self.parameters.binary_location
        )
        if not self._DOCK37_executor.is_available():
            raise DockingRunFailed("Cannot initialize DOCK 3.7 docker, as backend is not available - abort.")
        self._logger.log(f"Checked DOCK 3.7 backend availability (prefix_execution={self.parameters.prefix_execution}).",
                         _LE.DEBUG)

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetProp('dock37_score'))

    def add_molecules(self, molecules: list):
        """This method overrides the parent class, docker.py add_molecules method. This method appends prepared
        ligands to a list for subsequent docking. Note, DOCK 3.7 needs db2 format files.

        :param molecules: A list that is to contain all prepared ligands for subsequent docking
        :type molecules: list
        :raises NotImplementedError: Each backend must override the parent class, docker.py add_molecules method.
            Inability to do so or a bug causing incorrect implementation will raise a NotImplementedError
        """
        self.ligands = [mol.get_clone() for mol in molecules]
        self._docking_performed = False


    def _generate_temporary_input_output_files(self, start_indices, sublists):
        # in case singletons are handed over, wrap them in a list for "zipping" later
        if not isinstance(start_indices, list):
            start_indices = [start_indices]
        if not isinstance(sublists, list):
            sublists = [sublists]

        tmp_output_dirs = []
        db2_paths = []
        for start_index, sublist in zip(start_indices, sublists):
            cur_tmp_output_dir = tempfile.mkdtemp(dir='/tmp/qcxia02/', prefix='dock37')
            cur_tmp_docking_dir = f'{cur_tmp_output_dir}/docking'
            os.mkdir(cur_tmp_docking_dir)
            cur_tmp_input_db2_index = f'{cur_tmp_docking_dir}/split_database_index'
            with open(cur_tmp_input_db2_index, 'w') as sp_f:
                for ligand in sublist:
                    # write-out the temporary input file
                    if ligand.get_molecule() is None:
                        continue
                    db2_path = ligand.get_molecule()
                    sp_f.write(db2_path + '\n')
                    db2_paths.append(db2_path)
            tmp_output_dirs.append(cur_tmp_output_dir)
        return tmp_output_dirs, db2_paths


    def _dock(self, number_cores):

        self._initialize_executors()

        start_indices, sublists = self.get_sublists_for_docking(number_cores=number_cores)
        number_sublists = len(sublists) # 1 for 1 cpu cores (default)
        self._logger.log(f"Split ligands into {number_sublists} sublists for docking.", _LE.DEBUG)

        if not os.path.exists(self.parameters.dockfiles_indock_dir[0]):
            raise DockingRunFailed("Specified dockfiles and INDOCK directory to target (receptor) does not exist - abort.")

        sublists_submitted = 0
        slices_per_iteration = min(number_cores, number_sublists)
        while sublists_submitted < len(sublists):
            upper_bound_slice = min((sublists_submitted + slices_per_iteration), len(sublists))
            cur_slice_start_indices = start_indices[sublists_submitted:upper_bound_slice]
            cur_slice_sublists = sublists[sublists_submitted:upper_bound_slice]

            tmp_output_dirs, db2_paths = self._generate_temporary_input_output_files(cur_slice_start_indices, cur_slice_sublists)

            # run in parallel; wait for all subjobs to finish before proceeding
            processes = []
            for chunk_index in range(len(tmp_output_dirs)):
                p = multiprocessing.Process(target=self._dock_subjob, args=(tmp_output_dirs[chunk_index],))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()

            # add the number of input sublists rather than the output temporary folders to account for cases where
            # entire sublists failed to produce an input structure
            sublists_submitted += len(cur_slice_sublists)


            # parse the resulting sdf files
            for tmp_output_dir in tmp_output_dirs:
                # add conformations
                tmp_docked_poses = f'{tmp_output_dir}/poses.mol2'
                tmp_docked_unique_scores = f'{tmp_output_dir}/extract_all.sort.uniq.txt'
                tmp_docked_poses_unicon_sdf = f'{tmp_output_dir}/poses.mol2.sdf'
                if not os.path.isfile(tmp_docked_poses) or os.path.getsize(tmp_docked_poses) == 0:
                    continue

                os.system(f'{self.parameters.prefix_execution} && $UNICON_EXE -i {tmp_docked_poses} -o {tmp_docked_poses_unicon_sdf} &> /dev/null')

                if os.path.getsize(tmp_docked_poses_unicon_sdf) == 0:
                    self._logger.log(f"The size of {tmp_docked_poses_unicon_sdf} (converted by {tmp_docked_poses}) is 0.", _LE.WARNING)
                    continue

                for molecule in Chem.SDMolSupplier(tmp_docked_poses_unicon_sdf, removeHs=False):
                    if molecule is None:
                        continue

                    # extract the score from the Dock 3.7 output and update some tags
                    cur_identifier = molecule.GetProp('_Name').split()[0]
                    score = self._extract_score_from_Dock37Result(cur_identifier=cur_identifier, tmp_docked_unique_scores=tmp_docked_unique_scores)
                    molecule.SetProp("_Name", cur_identifier)
                    molecule.SetProp('dock37_score', score)
                    # molecule.ClearProp(_ROE.REMARK_TAG)

                    # add molecule to the appropriate ligand
                    for ligand in self.ligands:
                        if ligand.get_identifier() == cur_identifier:
                            ligand.add_conformer(molecule)
                            break

            # clean-up
            for path in tmp_output_dirs:
                shutil.rmtree(path)
            shutil.rmtree(os.path.dirname(os.path.dirname(db2_paths[0])))
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # the conformers are already sorted, but some tags are missing
        # -> <ligand_number>:<enumeration>:<conformer_number>
        for ligand in self.ligands:
            ligand.add_tags_to_conformers()

        # log any docking fails
        self._docking_fail_check()

        # generate docking results as dataframe
        result_parser = Dock37ResultParser(ligands=[ligand.get_clone() for ligand in self.ligands])
        self._df_results = result_parser.as_dataframe()

        # set docking flag
        self._docking_performed = True

    def _extract_score_from_Dock37Result(self, cur_identifier, tmp_docked_unique_scores) -> str:
        with open(tmp_docked_unique_scores, 'r') as score_f:
            score_lines = score_f.readlines()
        # result_line = [line for line in score_lines if cur_identifier in line][0]
        result_line = [line for line in score_lines if line and cur_identifier == line.split()[2]][0]
        parts = result_line.split()
        return parts[-1]

    def _dock_subjob(self, tmp_output_dir):

        # set up arguments list and execute
        # TODO: support "ensemble docking" - currently, only the first entry is used
        # search_space = self.parameters.search_space
        # arguments = [self.parameters.dockfiles_indock_dir[0], self.parameters.check_dir, f'&> {tmp_output_dir}/dock.log']
        arguments = [self.parameters.dockfiles_indock_dir[0], self.parameters.check_dir]


        execution_result = self._DOCK37_executor.execute(command='bash $DOCK37_SH',
                                                      arguments=arguments,
                                                      check=True, location=tmp_output_dir)
        self._logger.log(f"DOCK 3.7 output: {execution_result.stdout}.",
                         _LE.INFO)
        self._logger.log(f"DOCK 3.7 error: {execution_result.stderr}.",
                         _LE.INFO)
        if execution_result.returncode != 0:
            msg = f"Could not dock with DOCK 3.7, error message: {execution_result.stdout}."
            self._logger.log(msg, _LE.ERROR)
            raise DockingRunFailed(msg)
        
        tmp_docked_poses = f'{tmp_output_dir}/poses.mol2'
        self._delay4file_system(path=tmp_docked_poses)


    def write_docked_ligands(self, path, mode="all"):
        """This method overrides the parent class, docker.py write_docked_ligands method. This method writes docked
        ligands binding poses and conformers to a file. There is the option to output the best predicted binding pose
        per ligand, the best predicted binding pose per enumeration, or all the predicted binding poses

        :param path: Contains information on results output path
        :type path: string
        :param mode: Determines whether the output contains the best predicted binding pose per ligand, the best
            predicted binding pose per enumeration, or all the predicted binding poses
        :type mode: string, optional, default value is "all". Other possible values are "best_per_ligand" and
            "best_per_enumeration"
        :raises DockingRunFailed Error: This error is raised if the docking run has not been performed
        :raises OpenEye (OE) Fatal Error: This error is raised if the output file was unable to be created. Issues may
            be due to problems with the ligand structure
        :raises ValueError: This error is raised if the ligands are neither RDkit nor OpenEye readable
        """
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_RDKIT)
