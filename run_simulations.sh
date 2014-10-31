#!/bin/sh

TR=table3_raw.txt
rm -f $TR

NSIM=1000

for PAIRVAR in 0.1 0.2 0.5 1.0
do
    for SPECIESVAR in 0.1 0.2 0.5 1.0
    do
	for RESIDUALVAR in 0.1 0.2 0.5 1.0
	do
	    for EFFECTSIZE in 0.0 0.5
	    do
	    	Rscript evaluate_models.R $PAIRVAR $SPECIESVAR $RESIDUALVAR $EFFECTSIZE $NSIM 2> /dev/null | cut -f2- -d " " >> $TR
	    done
	done
    done
done
