#!/bin/sh
#
#PBS -N job
#PBS -l nodes=1:ppn=16
#PBS -l walltime=20:00:00
#PBS -j oe
source ~/.bash_profile
mkdir -p /scratch/ajasper
cd $PBS_O_WORKDIR
../../exe/nst.x < input > output
