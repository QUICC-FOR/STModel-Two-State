#!/bin/sh
DIR=~/STModel-Two-State
declare -a ScriptList=(ulmala quipri carova quecoc popgra jugnig ulmrub oxyarb ostvir franig)

for scr in ${ScriptList[@]}
do
	cd $DIR; qsub ./scr/3_anneal_queue/"$scr".sh
done
