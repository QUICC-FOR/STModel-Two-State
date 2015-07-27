#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N oxyarb

DIR=~/STModel-Two-State
SPECIES=23690-OXY-ARB
source $DIR/scr/3a_stm_model_selection.sh