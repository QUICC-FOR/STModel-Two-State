#!/bin/sh
#PBS -q default
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N pintae2

module load gcc/4.9.2
module load openmpi/1.8.3

# written for version 1.2
SRC=~/STModel-MCMC/bin


SPECIES=18037-PIN-TAE
DIR=~/STModel-Two-State/species/$SPECIES
cd $DIR; $SRC/stm2_mcmc -r res/mcmc2/resumeData.txt -t dat/mcmc_data.txt -i 80000 2>mcmc_log2.txt &

wait
