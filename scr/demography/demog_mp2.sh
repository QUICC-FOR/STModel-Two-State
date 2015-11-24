#!/bin/bash
#PBS -q qfat256
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -r n

module load R
cd /mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/mtalluto/STModel-Two-State

scr/demography/2_demog_compute.r