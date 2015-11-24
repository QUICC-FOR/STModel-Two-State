#!/bin/sh
#PBS -q default
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N demog
DIR=~/STModel-Two-State
cd $DIR; ./scr/demography/2_demog_compute.r

