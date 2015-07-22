#!/bin/sh
#PBS -q default
#PBS -l walltime=240:00:00
#PBS -l nodes=6:ppn=40
#PBS -r n
#PBS -N poptre.quealb.fraame23

module load gcc/4.9.2
module load openmpi/1.8.3

# written for version 1.2
SRC=~/STModel-MCMC/bin


SPECIES=195773-POP-TRE
DIR=~/STModel-Two-State/species/$SPECIES
cd $DIR; mkdir -p res/mcmc2 res/mcmc3
cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits2.txt -t dat/mcmc_data.txt -o res/mcmc2 -n 1 -b 0 -i 40000 -c 40 -l 5 -v 2 2>mcmc_log2.txt &
cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits3.txt -t dat/mcmc_data.txt -o res/mcmc3 -n 1 -b 0 -i 40000 -c 40 -l 5 -v 2 2>mcmc_log3.txt &

SPECIES=19290-QUE-ALB
DIR=~/STModel-Two-State/species/$SPECIES
cd $DIR; mkdir -p res/mcmc2 res/mcmc3
cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits2.txt -t dat/mcmc_data.txt -o res/mcmc2 -n 1 -b 0 -i 40000 -c 40 -l 5 -v 2 2>mcmc_log2.txt &
cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits3.txt -t dat/mcmc_data.txt -o res/mcmc3 -n 1 -b 0 -i 40000 -c 40 -l 5 -v 2 2>mcmc_log3.txt &

SPECIES=32931-FRA-AME
DIR=~/STModel-Two-State/species/$SPECIES
cd $DIR; mkdir -p res/mcmc2 res/mcmc3
cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits2.txt -t dat/mcmc_data.txt -o res/mcmc2 -n 1 -b 0 -i 40000 -c 40 -l 5 -v 2 2>mcmc_log2.txt &
cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits3.txt -t dat/mcmc_data.txt -o res/mcmc3 -n 1 -b 0 -i 40000 -c 40 -l 5 -v 2 2>mcmc_log3.txt &

wait
