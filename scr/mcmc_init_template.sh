

module load gcc/4.9.2
module load openmpi/1.8.3

SRC=~/STModel-MCMC/bin
DIR=~/STModel-Two-State/species/$SPECIES

# first step: run the MCMC initially for 8000 reps to get starting values
cd $DIR; $SRC/stm2_mcmc -p dat/mcmc_inits.txt -t dat/mcmc_data.txt -o res/mcmc1 -n 1 -b 0 -i 8000 -c 40 -l 5 -v 2 2>mcmc_log1.txt

# run an R script to generate new starting values
cd $DIR/../..; ./scr/generate_inits.r $SPECIES
