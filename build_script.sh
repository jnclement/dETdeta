#!/bin/bash

for ZCUT in 30; do #10 30
    for SCALE in 1; do #1 1.3 1.5
	for SUB in 0; do #0 0.018
	    for MINE in -10; do #0 0.005 0.018
		for SIM in 100; do
		    for DAT in 100; do
			for TAG in '_20231211_23696_etc_388p004'; do #'_20231113_nopileup_nzs' '_20231113_nopileup_wzs'; do
			    for TAG2 in '_20231211_etc'; do
				for COR in 1; do
				    root -b -q "build_hists.C(${SIM},${DAT},${ZCUT},${SCALE},${SUB},${MINE},\"${TAG}\",\"${TAG2}\",${COR})"
				done
			    done
			done
		    done
		done
	    done
	done
    done
done
