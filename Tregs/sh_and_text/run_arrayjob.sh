#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%A.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%A.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/ABM_CPP/Tregs
#SBATCH --export=NONE
###SBATCH --nodelist=eriskay,eorsa,harris,lewis,pabbay,sandray,vatersay,mingulay,canna
#SBATCH --exclude=eigg,tiree,raasay,taransay
#SBATCH --array=1-100
#SBATCH --cpus-per-task=5
 
while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

g++ -std=c++17 ${modelname}.cpp -o ${modelname} -larmadillo

echo "CODE Compiled"

time ./${modelname} $SLURM_ARRAY_TASK_ID

echo "END CODE"

