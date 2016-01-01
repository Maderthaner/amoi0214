
# useful web pages:
# psana-python: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual
# realtime plots: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-Real-timeOnlinePlotting/Monitoring
# mpi parallelization: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-MPIParallelization

# plotting commands (when run on the same node as the analysis):
# psplot OPAL &
# psplot XPROJ &

# plotting commands (when run on the a different node) '-s' means 'server':
# psplot -s daq-amo-mon02 OPAL &

# to get the analysis environment: source /reg/g/psdm/etc/ana_env.csh

# to run this on multiple nodes one needs to setup env properly on all of them by running the following (on line)
# /reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/mpirun -n 24 --host daq-amo-mon02,daq-amo-mon03,daq-amo-mon04 amoi0214-mpi.csh


# Standard PYTHON modules
print "IMPORTING STANDARD PYTHON MODULES...",
import matplotlib.pyplot as plt
import numpy as np
import math
import collections
import random
# print "DONE"

# LCLS psana to read data
# print "IMPORTING psana...",
import psana 
# print "DONE"

# For online plotting
from psmon import publish
from psmon.plots import XYPlot,Image

# custom algorithms
#from pypsalg import find_blobs
import find_blobs

# parallelization
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

try:
    opal_ped = np.load('ped.npy')
    if rank==0:
        print '*** Using pededestals from ped.npy'
except IOError:
    if rank==0:
        print '*** No pedestal subtraction will be done'

# custom functions
def moments(arr):
    if np.count_nonzero(arr) == 0:
        return 0,0
    nbins = len(arr)
    bins  = range(nbins)
    mean  = np.average(bins,weights=arr)
    var   = np.average((bins-mean)**2,weights=arr)
    sdev  = np.sqrt(var)
    return mean,sdev


if rank==0:
    publish.init()

# Set up detectors to read
opal_src = psana.Source("DetInfo(AmoEndstation.1:Opal1000.0)")
acq_src = psana.Source("DetInfo(AmoITOF.0:Acqiris.0)")
fee_src = psana.Source("BldInfo(FEEGasDetEnergy)")

# Set up event counters
eventCounter = 0
evtGood = 0
evtBad = 0
evtUsed = 0

# Buffers for histories
history_len = 6
history_len_long = 100
opal_hit_buff = collections.deque(maxlen=history_len)
opal_hit_avg_buff = collections.deque(maxlen=history_len_long)
opal_circ_buff = collections.deque(maxlen=history_len)
xproj_int_buff = collections.deque(maxlen=history_len)
xproj_int_avg_buff = collections.deque(maxlen=history_len)
moments_buff = collections.deque(maxlen=history_len_long)
#moments_avg_buff = collections.deque(maxlen=history_len_long)
#xhistogram_buff = collections.deque(maxlen=history_len)
hitxprojhist_buff = collections.deque(maxlen=history_len)
hitxprojjitter_buff = collections.deque(maxlen=history_len)
hitxprojseeded_buff = collections.deque(maxlen=history_len)

from Histogram import hist1d

hithist = hist1d(100,0.,1024.)
hitjitter = hist1d(100,0.,1024.)
hitseeded = hist1d(100,0.,1024.)

# ion yield array for sxrss scan
ion_yield = np.zeros(102) ##choose appropriate range

###### --- Online analysis
#ds = psana.DataSource('shmem=AMO.0:stop=no')
#for run in ds.runs():
#    for evt in run.events():
###### --- End of Online analysis region

###### --- Offline analysis
ds = psana.DataSource("exp=AMO/amoi0214:run=81:idx")
for run in ds.runs():
    times = run.times()
    mylength = len(times)/size
    mytimes= times[rank*mylength:(rank+1)*mylength]
    for i in range(mylength):
        evt = run.event(mytimes[i])
###### --- End of Offline analysis region

        #ebeam = evt.get(psana.Bld.BldDataEBeamV3,'BldInfo(Ebeam)')
        #ebeam = evt.get(psana.Bld.BldDataEBeamV6,'BldInfo(Ebeam)')
        #ebeam = evt.get(psana.Bld.EBeamL3E,'BldInfo(Ebeam)')
        gas_det = evt.get(psana.Bld.BldDataFEEGasDetEnergyV1,fee_src)
        #print str(gas_det.f_11_ENRC())
        #print str(gas_det.f_12_ENRC())
        print str(gas_det.f_21_ENRC())
        print str(gas_det.f_22_ENRC())