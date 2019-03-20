#!/bin/sh
#
#SBATCH --job-name=job
#SBATCH -N 1
#SBATCH --ntasks-per-node=36
#SBATCH -p knld
#SBATCH --time=20:00:00
#SBATCH -A KTP
#SBATCH -o job.%j.%N.out
source ~/.bash_profile
mkdir -p /scratch/ajasper
cd $SLURM_SUBMIT_DIR
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ajasper/lib/
export LD_LIBRARY_PATH
module load molpro
../../exe/nst.x m < input > output
