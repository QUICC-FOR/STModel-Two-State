#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N pintae2
SPECIES=18037-PIN-TAE

module load gcc/4.9.2
module load openmpi/1.8.3

# written for version 1.2
SRC=~/STModel-MCMC/bin
DIR=~/STModel-Two-State/species/$SPECIES

cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits2.txt -t dat/mcmc_data.txt -o res/mcmc2 -n 1 -b 0 -i 20000 -c 40 -l 5 -v 2 2>mcmc_log2.txt &

wait
