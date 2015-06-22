#!/bin/sh
DIR=~/STModel-Two-State
declare -a ScriptList=(nyssyl acerub lirtul quevel quenig pintae frapen pinban pinstr betall)

for scr in ${ScriptList[@]}
do
	cd $DIR; qsub ./scr/queue/"$scr".sh
done
