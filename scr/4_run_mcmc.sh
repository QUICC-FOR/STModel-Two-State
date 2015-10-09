#!/bin/bash
for scr in scr/mcmc/*; do
	qsub $scr
done