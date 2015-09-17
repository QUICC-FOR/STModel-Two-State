

module load gcc/4.9.2
module load openmpi/1.8.3

SRC=~/STModel-MCMC/bin
DIR=~/STModel-Two-State/species/$SPECIES

cd $DIR; $SRC/stm2_mcmc "$FLAG" -p dat/mcmc_inits"$N".txt -t dat/mcmc_data.txt -o res/mcmc"$N" -n 1 -b 0 -i 50000 -c 40 -l 5 -v 2 2>mcmc_log"$N".txt

# run an R script to generate new starting values
cd $DIR/../..; ./scr/generate_inits.r $SPECIES
