#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=30
#PBS -r n
#PBS -N lirtul

DIR=~/STModel-Two-State
SPECIES=18086-LIR-TUL
source $DIR/scr/3a_stm_model_selection.sh