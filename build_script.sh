#!/bin/bash

for ZCUT in 10 30; do
    for SCALE in 1 1.3; do
	for SUB in 0 0.018; do
	    for MINE in 0 0.005; do
		root -b -q "build_hists.C(${ZCUT},${SCALE},${SUB},${MINE})"
	    done
	done
    done
done
