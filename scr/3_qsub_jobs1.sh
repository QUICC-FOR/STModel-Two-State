#!/bin/sh
DIR=~/STModel-Two-State
declare -a ScriptList=(abibal acesac picgla picmar ulmame fraame quealb querub faggra poptre betpap)

for scr in ${ScriptList[@]}
do
	cd $DIR; qsub ./scr/queue/"$scr".sh
done
