#!/bin/sh

# this script launches a bunch of qsub jobs for the model selection/annealing for all species
# the output of this script is the entire res/anneal/ directory within each species dir

DIR=~/STModel-Two-State
#declare -a ScriptList=(ulmala quipri carova quecoc popgra jugnig ulmrub oxyarb ostvir franig)

for scr in `/bin/ls $DIR/scr/3_anneal_queue/`
do
	cd $DIR; qsub ./scr/3_anneal_queue/"$scr"
done
