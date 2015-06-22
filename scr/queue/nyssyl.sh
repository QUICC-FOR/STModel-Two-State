#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=30
#PBS -r n
#PBS -N nyssyl

DIR=~/STModel-Two-State
SPECIES=27821-NYS-SYL
source $DIR/scr/3a_stm_model_selection.sh