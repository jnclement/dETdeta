#!/bin/bash

for ZCUT in 30; do #10 30
    for SCALE in 1; do #1 1.3 1.5
	for SUB in 0; do #0 0.018
	    for MINE in -10; do #0 0.005 0.018
		for SIM in 1; do
		    for DAT in 1; do
			for TAG in '_20240226_data_run23696'; do
			    for TAG2 in '_20240226_nn_run23696'; do
				for RW in 0 1; do
				    for ZLO in `seq -15 13`; do
					ZUP=$(( $ZLO + 2 ))
					NAME="condor_${ZCUT}${SCALE}${SUB}${MINE}${DAT}${DAT}${TAG}${TAG2}${RW}${ZLO}${ZUP}.sh"
					echo "#!/bin/bash" > $NAME
					echo "source /opt/sphenix/core/bin/sphenix_setup.sh -n" >> $NAME
					echo "source /opt/sphenix/core/bin/setup_local.sh \"/sphenix/user/jocl/projects/testinstall\"" >> $NAME
					echo "export TESTINSTALL=\"/sphenix/user/jocl/projects/testinstall\"" >> $NAME
					echo "export HOME=/sphenix/u/jocl" >> $NAME
					echo "root -b -q \"build_hists.C(${DAT},${DAT},${ZCUT},${SCALE},${SUB},${MINE},\\\"${TAG}\\\",\\\"${TAG2}\\\",${RW},${ZLO},${ZUP})\"" >> $NAME
					chmod +x $NAME
					SUBN="condor_${ZCUT}${SCALE}${SUB}${MINE}${DAT}${DAT}${TAG}${TAG2}${RW}${ZLO}${ZUP}.sub"
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
	done
    done
done
