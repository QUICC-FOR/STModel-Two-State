#!/bin/bash
for dir in res/mcmc/*; do
	mkdir -p "$dir"/i0 "$dir"/ia "$dir"/ig	
done

for scr in scr/mcmc/int/*; do
	qsub $scr
done