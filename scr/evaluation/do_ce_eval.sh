#!/bin/sh
#PBS -q default
#PBS -l walltime=168:00:00
#PBS -l nodes=1:ppn=40
#PBS -r n
#PBS -N abibal-0
DIR=~/STModel-Two-State

cd $DIR; ./scr/ce_evaluation.r $SP $MOD 40 1000