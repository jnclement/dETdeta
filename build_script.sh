#!/bin/bash

for ZCUT in 30; do #10 30
    for SCALE in 1.5; do #1 1.3 1.5
	for SUB in 0; do #0 0.018
	    for MINE in 0 0.005 0.018; do #0 0.005 0.018
		root -b -q "build_hists.C+(${ZCUT},${SCALE},${SUB},${MINE})"
	    done
	done
    done
done
