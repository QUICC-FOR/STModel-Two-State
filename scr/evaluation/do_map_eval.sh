#!/bin/sh
#PBS -q default
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N mapeval

DIR=~/STModel-Two-State
MODELS=(0 a g i0 ia ig)

cd $DIR; ./scr/evaluation/map_eval.r $SP 0 40 1000
cd $DIR; ./scr/evaluation/map_eval.r $SP g 40 1000
cd $DIR; ./scr/evaluation/map_eval.r $SP i0 40 1000
cd $DIR; ./scr/evaluation/map_eval.r $SP ig 40 1000