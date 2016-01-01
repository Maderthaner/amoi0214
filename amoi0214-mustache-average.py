
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
# /reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/mpirun -n 60 --host psanacs061,psanacs062,psanacs060,psanacs059,psanacs058,psanacs057,psanacs056,psanacs055 amoi0214-spectra.csh

# to run on the batch farm, e.g. for run 82
#bsub -a mympi -n 2 -o junk.log -q psanaq python amoi0214-spectra.py 82

# Standard PYTHON modules
import numpy as np
import math
#import collections

# For offline plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# LCLS psana to read data
import psana

# custom algorithms
#from pypsalg import find_blobs
import find_blobs

# parallelization
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#rank = 0
#size = 1

try:
    opal_ped = np.load('ped.npy')
    if rank==0:
        print '*** Using pededestals from ped.npy'
except IOError:
    if rank==0:
        print '*** No pedestal subtraction will be done'

# For Histograms
from Histogram import hist1d
from pypsalg.Histogram import hist2d


# to pass commands to the python call line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("run", help="psana run number")
args = parser.parse_args()

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



# Set up detectors to read
if int(args.run) > 154:
    opal_src = psana.Source("DetInfo(AmoEndstation.0:Opal1000.0)") # For end of December runs
    print 'Using Opal Camera 0'
else:
    opal_src = psana.Source("DetInfo(AmoEndstation.1:Opal1000.0)") # For data of first runs
    print 'Using Opal Camera 1'
ebeam_src = psana.Source("BldInfo(EBeam)")

# Set up event counters
eventCounter = 0
evtGood = 0
evtBad = 0
evtUsed = 0

# Lists for mustacheplot
x = []
y = []

ebeamAvg = 0.0
hitCounterAvg = 0.0

if rank==0:
    print "DONE READING THE HEADER. CONTINUING ANALYSIS ON", size, "CORES FOR RUN", args.run,"."

###### --- Online analysis
#ds = psana.DataSource('shmem=AMO.0:stop=no')
#for run in ds.runs():
#    for evt in run.events():
###### --- End of Online analysis region

###### --- Offline analysis - initial fast call to get average data information
ds = psana.DataSource("exp=amoi0214:run=%s:idx" % args.run)
for run in ds.runs():
    if len(run.times()) > 5000:
        times = run.times()[0:5000]
    else:
        times = run.times()
    mylength = len(times)/size
    mytimes= times[rank*mylength:(rank+1)*mylength]
    for i in range(mylength):
        evt = run.event(mytimes[i])
###### --- End of Offline analysis region

        opal_raw = evt.get(psana.Camera.FrameV1,opal_src)
        ebeam = evt.get(psana.Bld.BldDataEBeamV6,ebeam_src)

        if opal_raw is None or ebeam is None:
            print 'Warning: Shot lacks data from Opal or eBeam.'
            continue
        else:
            evtUsed += 1

        ####### Copy detector data for further analysis #########
        opal = opal_raw.data16().astype(float)
        if 'opal_ped' in locals():
            opal -= opal_ped

        # threshold the image
        threshold = 700
        opal_thresh = (opal>threshold)*opal

        # find blobs
        c,w = find_blobs.find_blobs(opal,threshold)

        ebeamAvg += ebeam.ebeamL3Energy()
        hitCounterAvg += len(c)

        # Debugging
        if rank == 0 and evtUsed%10:
            print evtUsed*size

evtUsed_array = np.array([evtUsed])
if not 'evtUsed_sum_all' in locals():
    evtUsed_sum_all = np.zeros_like(evtUsed_array)
comm.Allreduce(evtUsed_array,evtUsed_sum_all)

ebeamAvg_array = np.array([ebeamAvg])
if not 'ebeamAvg_sum_all' in locals():
    ebeamAvg_sum_all = np.zeros_like(ebeamAvg_array)
comm.Allreduce(ebeamAvg_array,ebeamAvg_sum_all)

hitCounterAvg_array = np.array([hitCounterAvg])
if not 'hitCounterAvg_sum_all' in locals():
    hitCounterAvg_sum_all = np.zeros_like(hitCounterAvg_array)
comm.Allreduce(hitCounterAvg_array,hitCounterAvg_sum_all)

if rank == 0:
    print 'sum all:',evtUsed_sum_all[0],ebeamAvg_sum_all[0],hitCounterAvg_sum_all[0]

if rank==0 and evtUsed_sum_all[0] == 0:
    print 'Error: No good hits could be found in the initial 5000 hits for this run.'
    quit()

if rank==0 and ebeamAvg_sum_all[0] == 0:
    print 'Error: Average eBeam values could not be obtained'
    quit()

ebeamAvgNorm = (ebeamAvg_sum_all[0])/(evtUsed_sum_all[0])
hitCounterAvgNorm = (hitCounterAvg_sum_all[0])/(evtUsed_sum_all[0])

if hitCounterAvgNorm < 1:
    if rank == 0:
        print 'Warning: Average hitCounter value is less than one. Setting it artificaly to 1.'
    hitCounterAvgNorm = 1


ebeamLow = ebeamAvgNorm-20
ebeamHigh = ebeamAvgNorm+20
hitHigh = 2*hitCounterAvgNorm



if ebeamLow is None or ebeamHigh is None or hitHigh is None:
    if rank == 0:
        print 'Average values could not be obtained'
    quit()

if hitHigh < 100:
    hitHighSteps = np.ceil(1.5*np.ceil(hitHigh))
else:
    hitHighSteps = 100

if rank == 0:
    print 'avg. values',ebeamLow,ebeamHigh,hitHigh,hitHighSteps
    print 'Done with initial analysis to get average values.'

# initializing histograms
hithist = hist1d(100,0.,1023.)
hithistcorr = hist1d(100,0.,1023.)
mustachplot = hist2d(hitHighSteps,0,hitHigh,100,ebeamLow,ebeamHigh)
hithist1 = hist1d(100,0.,1023.)
hithist2 = hist1d(100,0.,1023.)
hithist3 = hist1d(100,0.,1023.)
hithist4 = hist1d(100,0.,1023.)
hithist5 = hist1d(100,0.,1023.)
hithist6 = hist1d(100,0.,1023.)
hithist7 = hist1d(100,0.,1023.)
hithist8 = hist1d(100,0.,1023.)
hithist9 = hist1d(100,0.,1023.)


###### --- Offline analysis
ds = psana.DataSource("exp=amoi0214:run=%s:idx" % args.run)
for run in ds.runs():
    times = run.times()
    mylength = len(times)/size
    mytimes= times[rank*mylength:(rank+1)*mylength]
    for i in range(mylength):
        evt = run.event(mytimes[i])
###### --- End of Offline analysis region

        opal_raw = evt.get(psana.Camera.FrameV1,opal_src)
        ebeam = evt.get(psana.Bld.BldDataEBeamV6,ebeam_src)

        #print ebeam.ebeamL3Energy()

        # Check all detectors are read in
        eventCounter += 1
        if opal_raw is None or ebeam is None:
            evtBad += 1
            #mustachplot.fill(0,3380)
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


        #################### DO BLOB FINDING #######################
        #### find blobs takes two variables. an Opal image (first) and then a threshold value (second)
        c,w = find_blobs.find_blobs(opal,threshold)

        ### Draw the blobs (careful, interupts the loop and plotting is 'slow' (Also, not MPI compatible...))
        #find_blobs.draw_blobs(opal,c,w) # draw the blobs in the opal picture


        #####find center of gravity. If there is a hit
        if len(c) != 0:
            x2,y2 = zip(*c)
            blob_centroid = np.sum(x2)/float(len(x2))
            #print blob_centroid
            shift = 512-int(blob_centroid)
        else:
            shift = 0
        
        ############## Spectra ################## 
        ######## Average Spectra
        #### Ebeam correlation plots 
        ### Row 1
        # Quadrant 1
        if ebeamLow < ebeam.ebeamL3Energy() < ((ebeamHigh-ebeamLow)/3)+ebeamLow and 0 < len(c) < hitHigh/3:
            for hit in c:
                hithist1.fill(hit[0])

        # Quadrant 2
        if ebeamLow < ebeam.ebeamL3Energy() < ((ebeamHigh-ebeamLow)/3)+ebeamLow and hitHigh/3 < len(c) < 2*hitHigh/3:
            for hit in c:
                hithist2.fill(hit[0])

        # Quadrant 3
        if ebeamLow < ebeam.ebeamL3Energy() < ((ebeamHigh-ebeamLow)/3)+ebeamLow and 2*hitHigh/3 < len(c) < 3*hitHigh/3:
            for hit in c:
                hithist3.fill(hit[0])

        ### Row 2
        # Quadrant 4
        if ((ebeamHigh-ebeamLow)/3)+ebeamLow < ebeam.ebeamL3Energy() < 2*((ebeamHigh-ebeamLow)/3)+ebeamLow and 0 < len(c) < hitHigh/3:
            for hit in c:
                hithist4.fill(hit[0])

        # Quadrant 5
        if ((ebeamHigh-ebeamLow)/3)+ebeamLow < ebeam.ebeamL3Energy() < 2*((ebeamHigh-ebeamLow)/3)+ebeamLow and hitHigh/3 < len(c) < 2*hitHigh/3:
            for hit in c:
                hithist5.fill(hit[0])

        # Quadrant 6
        if ((ebeamHigh-ebeamLow)/3)+ebeamLow < ebeam.ebeamL3Energy() < 2*((ebeamHigh-ebeamLow)/3)+ebeamLow and 2*hitHigh/3 < len(c) < 3*hitHigh/3:
            for hit in c:
                hithist6.fill(hit[0])

        ### Row 3
        # Quadrant 7
        if 2*((ebeamHigh-ebeamLow)/3)+ebeamLow < ebeam.ebeamL3Energy() < 3*((ebeamHigh-ebeamLow)/3)+ebeamLow and 0 < len(c) < hitHigh/3:
            for hit in c:
                hithist7.fill(hit[0])

        # Quadrant 8
        if 2*((ebeamHigh-ebeamLow)/3)+ebeamLow < ebeam.ebeamL3Energy() < 3*((ebeamHigh-ebeamLow)/3)+ebeamLow and hitHigh/3 < len(c) < 2*hitHigh/3:
            for hit in c:
                hithist8.fill(hit[0])
                
        # Quadrant 9
        if 2*((ebeamHigh-ebeamLow)/3)+ebeamLow < ebeam.ebeamL3Energy() < 3*((ebeamHigh-ebeamLow)/3)+ebeamLow and 2*hitHigh/3 < len(c) < 3*hitHigh/3:
            for hit in c:
                hithist9.fill(hit[0])

        # Mustacheplot
        mustachplot.fill(len(c),ebeam.ebeamL3Energy())

        x.append(ebeam.ebeamL3Energy())
        y.append(len(c))

        # Debugging
        if rank == 0 and evtGood%10:
            print evtGood*size

############### create zero arrays and dump for master

xall = comm.allgather(x)
yall = comm.allgather(y)

if not 'mustachplot_sum_all' in locals():
    mustachplot_sum_all = np.empty_like(mustachplot.data)
comm.Allreduce(mustachplot.data,mustachplot_sum_all)

if not 'hithist1_sum_all' in locals():
    hithist1_sum_all = np.empty_like(hithist1.data)
comm.Allreduce(hithist1.data,hithist1_sum_all)

if not 'hithist2_sum_all' in locals():
    hithist2_sum_all = np.empty_like(hithist2.data)
comm.Allreduce(hithist2.data,hithist2_sum_all)

if not 'hithist2_sum_all' in locals():
    hithist2_sum_all = np.empty_like(hithist2.data)
comm.Allreduce(hithist2.data,hithist2_sum_all)

if not 'hithist3_sum_all' in locals():
    hithist3_sum_all = np.empty_like(hithist3.data)
comm.Allreduce(hithist3.data,hithist3_sum_all)

if not 'hithist4_sum_all' in locals():
    hithist4_sum_all = np.empty_like(hithist4.data)
comm.Allreduce(hithist4.data,hithist4_sum_all)

if not 'hithist5_sum_all' in locals():
    hithist5_sum_all = np.empty_like(hithist5.data)
comm.Allreduce(hithist5.data,hithist5_sum_all)

if not 'hithist6_sum_all' in locals():
    hithist6_sum_all = np.empty_like(hithist6.data)
comm.Allreduce(hithist6.data,hithist6_sum_all)

if not 'hithist7_sum_all' in locals():
    hithist7_sum_all = np.empty_like(hithist7.data)
comm.Allreduce(hithist7.data,hithist7_sum_all)

if not 'hithist8_sum_all' in locals():
    hithist8_sum_all = np.empty_like(hithist8.data)
comm.Allreduce(hithist8.data,hithist8_sum_all)

if not 'hithist9_sum_all' in locals():
    hithist9_sum_all = np.empty_like(hithist9.data)
comm.Allreduce(hithist9.data,hithist9_sum_all)

evtGood_array = np.array([evtGood])
if not 'evtGood_sum_all' in locals():
    evtGood_sum_all = np.zeros_like(evtGood_array)
comm.Allreduce(evtGood_array,evtGood_sum_all)


if rank==0:
    ###### Output (debugging and info)
    print evtGood_sum_all[0], 'total events processed. Saving now...' # Real number of good events
    
    ##### Saving spectra to textfiles #####
    #names
    mustacheplot_name = 'mustacheplot/run'+str(args.run).zfill(4)+'_MustachePlot'
    quadrant_name = 'mustacheplot/run'+str(args.run).zfill(4)+'_BinnedSpectraQuadrant'
    
    #saveing files
    np.savetxt(mustacheplot_name+'_xVals',xall)
    np.savetxt(mustacheplot_name+'_yVals',yall)
    np.savetxt(mustacheplot_name+'_totalHitCount',evtGood_sum_all)
    np.savetxt(mustacheplot_name,mustachplot_sum_all)
    np.savetxt(quadrant_name+'_1',hithist1_sum_all)
    np.savetxt(quadrant_name+'_2',hithist2_sum_all)
    np.savetxt(quadrant_name+'_3',hithist3_sum_all)
    np.savetxt(quadrant_name+'_4',hithist4_sum_all)
    np.savetxt(quadrant_name+'_5',hithist5_sum_all)
    np.savetxt(quadrant_name+'_6',hithist6_sum_all)
    np.savetxt(quadrant_name+'_7',hithist7_sum_all)
    np.savetxt(quadrant_name+'_8',hithist8_sum_all)
    np.savetxt(quadrant_name+'_9',hithist9_sum_all)

    print 'Files have been written into the working_directory/mustacheplot/.'
    
    ############### Exquisit plotting##################
    with PdfPages(quadrant_name+'.pdf') as pdf:
        fig = plt.figure()
        #plt.subplot(2,1,1)
        #plt.plot(range(0,len(mustachplot_sum_all)),mustachplot_sum_all/evtGood_sum_all[0])
        #plt.xlabel('Ebeam L3 energy in eV')
        #plt.ylabel('Intensity in arb. units')
        # plotting Row 1
        plt.subplot(3,1,1)
        plt.plot(range(0,len(hithist1_sum_all)),hithist1_sum_all)
        plt.subplot(3,1,1)
        plt.plot(range(0,len(hithist2_sum_all)),hithist2_sum_all)
        plt.subplot(3,1,1)
        plt.plot(range(0,len(hithist3_sum_all)),hithist3_sum_all)
        plt.xlabel('Bins')
        plt.ylabel('Intensity in arb. units')
        # plotting row 2
        plt.subplot(3,1,2)
        plt.plot(range(0,len(hithist4_sum_all)),hithist4_sum_all)
        plt.subplot(3,1,2)
        plt.plot(range(0,len(hithist5_sum_all)),hithist5_sum_all)
        plt.subplot(3,1,2)
        plt.plot(range(0,len(hithist6_sum_all)),hithist6_sum_all)
        plt.xlabel('Bins')
        plt.ylabel('Intensity in arb. units')
        # plotting row 3
        plt.subplot(3,1,3)
        plt.plot(range(0,len(hithist7_sum_all)),hithist7_sum_all)
        plt.subplot(3,1,3)
        plt.plot(range(0,len(hithist8_sum_all)),hithist8_sum_all)
        plt.subplot(3,1,3)
        plt.plot(range(0,len(hithist9_sum_all)),hithist9_sum_all)
        plt.xlabel('Bins')
        plt.ylabel('Intensity in arb. units')
        pdf.savefig(fig)
        print 'Plot succesful written into working_directory/mustacheplot/.'

    print 'Analysis completed'

                