#!/bin/sh
#PBS -q default
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N mapeval
DIR=~/STModel-Two-State
cd $DIR; ./scr/evaluation/map_eval.r $SP $MOD 40 1000
