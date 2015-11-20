#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=10:ppn=40
#PBS -r n
#PBS -N resp_curves
SPECIES=(18032-ABI-BAL 19049-ULM-AME 19290-QUE-ALB 19408-QUE-RUB 19481-BET-ALL 19489-BET-PAP 28728-ACE-RUB 28731-ACE-SAC 183302-PIC-MAR 195773-POP-TRE)
DIR=~/STModel-Two-State
for SP in "${SPECIES[@]}"
do
    cd $DIR; ./scr/post_processing/b_compute_resp_curves.r $SP 40 &
done
wait
