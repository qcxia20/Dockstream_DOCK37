#!/bin/bash


insmi=$1
max_conf=$2
workingpath=$3
method=$4
checkstereo=$5
useff=$6
sampletp=$7
reseth=$8
rotateh=$9
keep_max_conf=${10}
mergeiso=${11}

# source $CONDA_ROOT/etc/profile.d/conda.sh
source /pubhome/qcxia02/mambaforge/etc/profile.d/conda.sh
conda activate DockStreamCommunity

if [ $sampletp == 'True' ]; then
    $UNICON_EXE -i $insmi -p single -t single -o $insmi.uni.smi &> /dev/null
    mv $insmi.uni.smi $insmi
fi

export OMP_NUM_THREADS=1
command="build_ligand -i $insmi -n $max_conf --workingpath $workingpath --outputpath $workingpath --method $method"
if [ $checkstereo == 'True' ];
then
    command+=" -c"
fi

if [ $useff == 'True' ];
then
    command+=" -f"
fi

if [ $reseth == 'True' ];
then
    command+=" --reseth"
fi

if [ $rotateh == 'True' ];
then
    command+=" --rotateh"
fi

if [ $keep_max_conf == 'True' ];
then
    command+=" --keep_max_conf"
fi

if [ $mergeiso == 'True' ];
then
    command+=" --mergeiso"
fi
# echo $command

# echo $command
$command &> db2_${method}_$(basename $insmi .smi).log
