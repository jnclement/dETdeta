#!/bin/bash
NO=0
for i in {0..200}; do
    root -b -q $1 2>&1 > output${i}.txt
    if [[ $i > 0 ]]
    then
	OTH=$(( $i - 1 ))
	diff output${i}.txt output${OTH}.txt
	RET=$?
	if [[ $RET == "0" ]]
	then
	    echo "Same"
	else
	    echo "No!"
	    NO=$(( $NO + 1 ))
	fi
    fi
    if [[ $NO > 10 ]]
    then
	echo ">5 unsame found."
	break
    fi
done
