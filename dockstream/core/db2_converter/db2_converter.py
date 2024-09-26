import os
import glob
import tempfile
from copy import deepcopy
import multiprocessing
from pathlib import Path

from typing import Optional, List
from pydantic import PrivateAttr, BaseModel
from rdkit import Chem
from typing_extensions import Literal

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed
from dockstream.utils.execute_external.db2_converter import DB2_converter_Executor

from dockstream.core.ligand_preparator import LigandPreparator, _LE
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.smiles import to_mol
from dockstream.core.ligand.ligand import Ligand

_LP = RDkitLigandPreparationEnum()

# This class is inspired / based on code written by
# Peter Schmidtke (https://github.com/Discngine/rdkit_tethered_minimization), which is accessible under the MIT license
# and the blogspot entry called "more-on-constrained-embedding" from the RDkit guys.


class DB2_converterParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    max_conf: Optional[int] = 1000
    db2_method: Optional[str] = 'conformator'
    checkstereo: Optional[bool] = False
    useff: Optional[bool] = False
    reseth: Optional[bool] = False
    rotateh: Optional[bool] = False
    keep_max_conf: Optional[bool] = False
    sampletp: Optional[bool] = False
    mergeiso: Optional[bool] = False


class DB2_converter(LigandPreparator, BaseModel):
    """Class that deals with all the preparatory steps needed before actual docking using "rDock" can commence."""

    type: Literal["db2_converter"] = "db2_converter"
    parameters: DB2_converterParameters = DB2_converterParameters()

    _db2_converter_executor: DB2_converter_Executor = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)
        self._db2_converter_executor = DB2_converter_Executor(prefix_execution=self.parameters.prefix_execution,
                                                              binary_location=self.parameters.binary_location)
        if not self._db2_converter_executor.is_available():
            raise LigandPreparationFailed("Cannot initialize DB2_converter backend - abort.")

    def _load_references(self):
        references = []
        ref_format = self.align.reference_format.upper()
        for path in self.align.reference_paths:
            if ref_format == _LP.ALIGN_REFERENCE_FORMAT_PDB:
                ref_mol = Chem.MolFromPDBFile(path, sanitize=True)
                ref_mol.SetProp("_Name", os.path.basename(path))
                references.append(ref_mol)
            elif ref_format == _LP.ALIGN_REFERENCE_FORMAT_SDF:
                mol_supplier = Chem.SDMolSupplier(path)
                for mol in mol_supplier:
                    references.append(mol)
            else:
                raise IOError("Specified format not supported.")
        if len(references) == 0:
            raise LigandPreparationFailed("No reference molecules could be loaded with path(s) specified.")
        self._references = references
        self._logger.log(f"Stored {len(references)} reference molecules.", _LE.DEBUG)

    def _smiles_to_molecules(self, ligands: List[Ligand]) -> List[Ligand]:
        for lig in ligands:
            mol = to_mol(lig.get_smile())
            lig.set_molecule(mol)
            lig.set_mol_type(_LP.TYPE_RDKIT)
        return ligands

    def generate3Dcoordinates(self):
        """Method to generate 3D coordinates, in case the molecules have been built from SMILES."""

        for lig in self.ligands:
            lig.set_molecule(None)
            lig.set_mol_type(None)
        ligand_list = self._smiles_to_molecules(deepcopy(self.ligands))

        failed = 0
        succeeded = 0

        tmppath = f"/tmp/qcxia02/"
        # print("tmppath", tmppath)
        if not os.path.exists(tmppath):
            os.mkdir(tmppath)
        tmp_dir = tempfile.mkdtemp(dir=tmppath)
        all_smi_path = f'{tmp_dir}/all.smi'
        with open(all_smi_path, 'w') as f:
            for lig in self.ligands:
                f.write(lig.get_smile() + " " + lig.get_identifier() + "\n")

        # split all.smi into cpu cores
        Ncores = multiprocessing.cpu_count()
        Nsmis = len([ line for line in Path(all_smi_path).read_text().split("\n") if line ])
        Nsplit = Nsmis // Ncores + 1
        os.system(f"split -l {Nsplit} {all_smi_path} {tmp_dir}/split. --additional-suffix=.smi")
        # os.system
        db2_converter_command=f"parallel -j {Ncores} -k 'timeout 300 bash $db2_converter_SH {{}} {self.parameters.max_conf} {tmp_dir}/input {self.parameters.db2_method} {self.parameters.checkstereo} {self.parameters.useff} {self.parameters.sampletp} {self.parameters.reseth} {self.parameters.rotateh} {self.parameters.keep_max_conf} {self.parameters.mergeiso}'"
        # print(db2_converter_command)
        os.system(f"{self.parameters.prefix_execution}; cwd=`pwd`; cd {tmp_dir}; ls split.*.smi | {db2_converter_command}; cd $cwd") # each 5 mins at most


        for idx, lig_obj in enumerate(ligand_list):
            if not os.path.exists(f'{tmp_dir}/input/{lig_obj.get_identifier()}.db2.gz'):
                if glob.glob(f'{tmp_dir}/input/{lig_obj.get_identifier()}.*.db2.gz'):
                    os.system(f"cat {tmp_dir}/input/{lig_obj.get_identifier()}.*.db2.gz > {tmp_dir}/input/{lig_obj.get_identifier()}.db2.gz")
                    os.system(f"rm {tmp_dir}/input/{lig_obj.get_identifier()}.*.db2.gz")
            if not os.path.exists(f'{tmp_dir}/input/{lig_obj.get_identifier()}.db2.gz'):
                self._logger.log(f"Could not embed molecule number {lig_obj.get_ligand_number()} (smile: {lig_obj.get_smile()}) - no 3D coordinates generated.",
                                    _LE.DEBUG)
                failed += 1
                continue

            self.ligands[idx] = Ligand(smile=lig_obj.get_smile(),
                                       original_smile=lig_obj.get_original_smile(),
                                       ligand_number=lig_obj.get_ligand_number(),
                                       enumeration=lig_obj.get_enumeration(),
                                       molecule=f'{tmp_dir}/input/{lig_obj.get_identifier()}.db2.gz',
                                       mol_type=lig_obj.get_mol_type(),
                                       name=lig_obj.get_name())
            succeeded += 1

        if failed > 0:
            self._logger.log(f"Of {len(self.ligands)}, {failed} could not be embedded.",
                             _LE.WARNING)
        self._logger.log(f"In total, {succeeded} ligands were successfully embedded (db2).", _LE.DEBUG)

    def align_ligands(self):
        self.ligands = self._align_ligands_with_RDkit_preparator(self.ligands)
