#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N ulmame3
SPECIES=19049-ULM-AME

module load gcc/4.9.2
module load openmpi/1.8.3

# written for version 1.2
SRC=~/STModel-MCMC/bin
DIR=~/STModel-Two-State/species/$SPECIES

cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits3.txt -t dat/mcmc_data.txt -o res/mcmc3 -n 1 -b 0 -i 20000 -c 40 -l 5 -v 2 2>mcmc_log3.txt &

wait