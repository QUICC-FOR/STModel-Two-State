#!/bin/sh
DIR=~/STModel-Two-State
declare -a ScriptList=(liqsty caralb queste cargla quefal tsucan junvir pinech thuocc picrub)

for scr in ${ScriptList[@]}
do
	cd $DIR; qsub ./scr/3_anneal_queue/"$scr".sh
done
