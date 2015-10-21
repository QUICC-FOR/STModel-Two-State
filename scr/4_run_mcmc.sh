#!/bin/bash
for scr in scr/mcmc/*; do
	qsub $scr
done

# run scripts for intercept only models
for scr in scr/mcmc/int/*; do
	qsub $scr
done