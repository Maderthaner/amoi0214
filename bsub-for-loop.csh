#!/bin/csh
foreach x (`seq 158 1 208`)
	bsub -a mympi -n 8 -o run_$x.log -q psanaq -R "span[ptile=1]" python amoi0214-spectra.py $x
	sleep 1
end
