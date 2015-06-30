#!/bin/sh
#PBS -q default
#PBS -l walltime=120:00:00
#PBS -l nodes=2:ppn=40
#PBS -r n
#PBS -N abibal.acesac.1

module load gcc/4.9.2
module load openmpi/1.8.3

# written for version 1.2
SRC=~/STModel-MCMC/bin


SPECIES=18032-ABI-BAL
DIR=~/STModel-Two-State/species/$SPECIES
cd $DIR; mkdir -p res/mcmc1
cd $DIR; $SRC/stm2_mcmc -r res/mcmc1/resumeData.txt -t dat/mcmc_data.txt -i 20000 2>mcmc_log1.txt &

SPECIES=28731-ACE-SAC
DIR=~/STModel-Two-State/species/$SPECIES
cd $DIR; mkdir -p res/mcmc1
cd $DIR; $SRC/stm2_mcmc -r res/mcmc1/resumeData.txt -t dat/mcmc_data.txt -i 20000 2>mcmc_log1.txt &

wait
