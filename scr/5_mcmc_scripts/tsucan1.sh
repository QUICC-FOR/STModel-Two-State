#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N tsucan1

module load gcc/4.9.2
module load openmpi/1.8.3

# written for version 1.2
SRC=~/STModel-MCMC/bin


SPECIES=183397-TSU-CAN
DIR=~/STModel-Two-State/species/$SPECIES
cd $DIR; $SRC/stm2_mcmc -r res/mcmc1/resumeData.txt -t dat/mcmc_data.txt -i 10000 2>mcmc_log1.txt &

wait
