#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N jugnig

DIR=~/STModel-Two-State
SPECIES=19254-JUG-NIG
source $DIR/scr/stm_model_selection_4.sh