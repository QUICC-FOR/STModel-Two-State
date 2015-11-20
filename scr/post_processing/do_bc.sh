#!/bin/sh
#PBS -q default
#PBS -l walltime=24:00:00
#PBS -l nodes=2:ppn=40
#PBS -r n
#PBS -N rc_grid
DIR=~/STModel-Two-State
cd $DIR; ./scr/post_processing/b_compute_resp_curves.r $SP 40 &
cd $DIR; ./scr/post_processing/c_compute_grid_preds.r $SP 40 &
wait
