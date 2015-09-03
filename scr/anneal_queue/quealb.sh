#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N quealb

DIR=~/STModel-Two-State
SPECIES=19290-QUE-ALB
source $DIR/scr/stm_model_selection_4.sh