#!/bin/bash
SPECIES=(18032-ABI-BAL 19049-ULM-AME 19290-QUE-ALB 19408-QUE-RUB 19481-BET-ALL 19489-BET-PAP 28728-ACE-RUB 28731-ACE-SAC 183302-PIC-MAR 195773-POP-TRE)
MODELS=(0 a g i0 ia ig)
DIR=~/STModel-Two-State
for sp in "${SPECIES[@]}"
do
	for mod in "${MODELS[@]}"
	do
 		qsub -vSP=$sp,MOD=$mod $DIR/scr/do_ce_eval.sh
	done
done