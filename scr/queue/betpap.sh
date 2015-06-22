#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=30
#PBS -r n
#PBS -N betpap

SPECIES=19489-BET-PAP
source scr/3a_stm_model_selection.sh