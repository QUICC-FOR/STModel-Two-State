#!/bin/sh
#PBS -q qwork
#PBS -l walltime=40:00:00
#PBS -l nodes=1:ppn=1
#PBS -r n

# note that entire nodes are always allocated, so
# for qwork nodes=1:ppn=1 gives 24 cores
# for qfat256 you get 48 cores
# qfbb has 24 per node, but the minimum request is 12 nodes

module load gcc/4.8.2
module load gsl64/1.16

BIN=~/STModel-MCMC/bin/stm2_1.5.1
DIR=~/STModel-Two-State

sp="505490-THU-OCC"
ni=100000

$BIN -r $DIR/res/mcmc/"$sp"/0_ch2/resumeData.txt -t $DIR/dat/mcmc/"$sp"_mcmcDat.txt -i "$ni" 2>$DIR/log/mcmc_log_"$sp"_0_ch2.txt &
	
wait
