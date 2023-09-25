#!/bin/bash

for ZCUT in 10 30; do #10 30
    for SCALE in 1 1.1 1.3 1.5; do #1 1.3 1.5
	for SUB in 0 0.01; do #0 0.018
	    for MINE in 0 0.05 0.01 0.018; do #0 0.005 0.018
		for SIM in 1; do
		    for DAT in 1; do
			for TAG in '_otc' '_ntc'; do
			    NAME="condor_${ZCUT}${SCALE}${SUB}${MINE}${SIM}${DAT}${TAG}.sh"
			    echo "#!/bin/bash" > $NAME
			    echo "source /opt/sphenix/core/bin/sphenix_setup.sh -n" >> $NAME
			    echo "source /opt/sphenix/core/bin/setup_local.sh \"/sphenix/user/jocl/projects/testinstall\"" >> $NAME
			    echo "export TESTINSTALL=\"/sphenix/user/jocl/projects/testinstall\"" >> $NAME
			    echo "export HOME=/sphenix/u/jocl" >> $NAME
			    echo "root -b -q \"build_hists.C(${SIM},${DAT},${ZCUT},${SCALE},${SUB},${MINE},\\\"${TAG}\\\")\"" >> $NAME
			    chmod +x $NAME
			    SUBN="condor_${ZCUT}${SCALE}${SUB}${MINE}${SIM}${DAT}${TAG}.sub"
			    echo "executable = ${NAME}" > $SUBN
			    echo "arguments =" >> $SUBN
			    echo "output = run/output/out/output${SUBN}.out" >> $SUBN
			    echo "should_transfer_files   = IF_NEEDED" >> $SUBN
			    echo "when_to_transfer_output = ON_EXIT" >> $SUBN
			    echo "error                   = run/output/err/error${SUBN}.err" >> $SUBN
			    echo "log = /tmp/jocl_${SUBN}.log" >> $SUBN
			    echo "queue 1" >> $SUBN
			    condor_submit $SUBN
			done
		    done
		done
	    done
	done
    done
done
