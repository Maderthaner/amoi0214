# useful web pages:
# psana-python: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual
# realtime plots: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-Real-timeOnlinePlotting/Monitoring
# mpi parallelization: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-MPIParallelization

# plotting commands (when run on the same node as the analysis):
# psplot OPAL &
# psplot XPROJ &

# plotting commands (when run on the a different node) '-s' means 'server':
# psplot -s daq-amo-mon03 OPAL &

# to get the analysis environment: source /reg/g/psdm/etc/ana_env.csh

# to run this on multiple nodes one needs to setup env properly on all of them by running the following (on line)
# /reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/mpirun -n 16 --host daq-amo-mon03,daq-amo-mon04 amoi0214-mpi.csh


# Standard PYTHON modules
print "STARTING...",
import matplotlib.pyplot as plt
import numpy as np
import math
import collections
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
opal_src = psana.Source("DetInfo(AmoEndstation.1:Opal1000.0)") # For data of first runs
#opal_src = psana.Source("DetInfo(AmoEndstation.0:Opal1000.0)") # End of December data

# Set up event counters
eventCounter = 0
evtGood = 0
evtBad = 0
evtUsed = 0

# Buffers for histories
history_len = 6
history_len_long = 50
opal_hit_buff = collections.deque(maxlen=history_len)
opal_hit_avg_buff = collections.deque(maxlen=history_len_long)
opal_circ_buff = collections.deque(maxlen=history_len)
xproj_max_buff = collections.deque(maxlen=history_len)
#xproj_int_avg_buff = collections.deque(maxlen=history_len)
moments_buff = collections.deque(maxlen=history_len_long)
#moments_avg_buff = collections.deque(maxlen=history_len_long)
#xhistogram_buff = collections.deque(maxlen=history_len)
hitxprojhist_buff = collections.deque(maxlen=history_len)

from Histogram import hist1d

hithist = hist1d(100,0.,1024.)
jitterhist = hist1d(100,-512.,512.)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("run", help="psana run number")
args = parser.parse_args()



if rank==0:
    publish.init()
    print "DONE READING THE HEADER. CONTINUING ANALYSIS ON", size, "CORES FOR RUN", args.run,"."



###### --- Online analysis
# ds = psana.DataSource('shmem=AMO.0:stop=no')
# for run in ds.runs():
#     for evt in run.events():
###### --- End of Online analysis region

###### --- Offline analysis
ds = psana.DataSource("exp=AMO/amoi0214:run=%s:idx" % args.run)
for run in ds.runs():
    times = run.times()
    mylength = len(times)/size
    mytimes= times[rank*mylength:(rank+1)*mylength]
    for i in range(mylength):
        evt = run.event(mytimes[i])
###### --- End of Offline analysis region

        opal_raw = evt.get(psana.Camera.FrameV1,opal_src)

        # Check all detectors are read in
        eventCounter += 1
        if (opal_raw is None) :
            evtBad += 1
            continue
        else:
            evtGood += 1

        ####### Copy detector data for further analysis #########
        opal = opal_raw.data16().astype(float)
        if 'opal_ped' in locals():
            opal -= opal_ped
            

        # threshold the image
        threshold = 700
        opal_thresh = (opal>threshold)*opal

        # do two projections
        opal_thresh_xproj = np.sum(opal_thresh,axis=1)

        # sum up the projected image in bin range 100 to 200
        #integrated_projx = np.sum(opal_thresh_xproj[100:200])

        # do blob finding
        # find blobs takes two variables. an Opal image (first) and then a threshold value (second)
        c,w = find_blobs.find_blobs(opal,threshold)

        #find center of gravity. If there is a hit
        if len(c) != 0:
            x2,y2 = zip(*c)
            blob_centroid = np.sum(x2)/float(len(x2))
            shift = 512-int(blob_centroid)
        else:
            shift = 0
            
        # Hit Histogram
        hithist.data = np.zeros(hithist.nbinx)
        for hit in c:
            hithist.fill(hit[0]+shift)

        # Jitter histogram
        jitterhist.fill(shift)


    if not 'jitterhist_sum_all' in locals():
        jitterhist_sum_all = np.empty_like(jitterhist.data)
    comm.Reduce(jitterhist.data,jitterhist_sum_all)

    if rank==0:

        spectra_name = 'jitterdata/run'+str(args.run).zfill(4)+'_jitterhist'
        np.savetxt(spectra_name,jitterhist_sum_all)

        print 'File has been written into the working_directory/jitterdata/.'
