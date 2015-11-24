#!/bin/sh
#PBS -q qfbb
#PBS -l walltime=8:00:00
#PBS -l nodes=12:ppn=1
#PBS -r n

# note that entire nodes are always allocated, so
# for qwork nodes=1:ppn=1 gives 24 cores
# for qfat256 you get 48 cores
# qfbb has 24 per node, but the minimum request is 12 nodes

module load gcc/4.8.2
module load gsl64/1.16

BIN=~/STModel-MCMC/bin/stm2_1.5.1
DIR=~/STModel-Two-State

SPLIST=(18034-PIC-RUB 183295-PIC-GLA 505490-THU-OCC 183319-PIN-BAN)

for sp in "${SPLIST[@]}"
do
	mkdir -p $DIR/res/mcmc/"$sp"/0i_ch1
	$BIN -d -p $DIR/dat/mcmc/"$sp"_mcmcInit_int1.txt -t $DIR/dat/mcmc/"$sp"_mcmcDat.txt -i 20000 -c 24 -l 5 -v 2 -o $DIR/res/mcmc/"$sp"/0i_ch1  2>$DIR/log/mcmc_log_"$sp"_0i_ch1.txt &
	
	mkdir -p $DIR/res/mcmc/"$sp"/0i_ch2
	$BIN -d -p $DIR/dat/mcmc/"$sp"_mcmcInit_int2.txt -t $DIR/dat/mcmc/"$sp"_mcmcDat.txt -i 20000 -c 24 -l 5 -v 2 -o $DIR/res/mcmc/"$sp"/0i_ch2  2>$DIR/log/mcmc_log_"$sp"_0i_ch2.txt &

	mkdir -p $DIR/res/mcmc/"$sp"/0i_ch3
	$BIN -d -p $DIR/dat/mcmc/"$sp"_mcmcInit_int3.txt -t $DIR/dat/mcmc/"$sp"_mcmcDat.txt -i 20000 -c 24 -l 5 -v 2 -o $DIR/res/mcmc/"$sp"/0i_ch3  2>$DIR/log/mcmc_log_"$sp"_0i_ch3.txt &
	
done
wait
