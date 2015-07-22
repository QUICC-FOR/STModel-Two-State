#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N pintae

DIR=~/STModel-Two-State
SPECIES=18037-PIN-TAE
source $DIR/scr/3a_stm_model_selection.sh