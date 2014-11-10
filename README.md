### useful web pages:

psana-python: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual

realtime plots: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-Real-timeOnlinePlotting/Monitoring

mpi parallelization: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-MPIParallelization

### plotting commands (when run on the same node as the analysis):

psplot OPAL &
psplot XPROJ &

### plotting commands (when run on the a different node) '-s' means 'server':
psplot -s daq-amo-mon03 OPAL &

### to get the analysis environment: source /reg/g/psdm/etc/ana_env.csh

to run this on multiple nodes one needs to setup env properly on all of them by running the following (on line)
/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/mpirun -n 16 --host daq-amo-mon03,daq-amo-mon04 amoi0214-mpi.csh
