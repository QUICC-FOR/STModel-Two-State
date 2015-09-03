#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N queste

DIR=~/STModel-Two-State
SPECIES=19422-QUE-STE
source $DIR/scr/stm_model_selection_4.sh