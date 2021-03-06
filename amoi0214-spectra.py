
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
import matplotlib.pyplot as plt
import numpy as np
import math
#import collections
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
opal_src = psana.Source("DetInfo(AmoEndstation.1:Opal1000.0)")

# Set up event counters
eventCounter = 0
evtGood = 0
evtBad = 0
evtUsed = 0

# initializing histograms
hithist = hist1d(100,0.,1023.)
hithistcorr = hist1d(100,0.,1023.)


if rank==0:
    publish.init()
    print "DONE READING THE HEADER. CONTINUING ANALYSIS ON", size, "CORES FOR RUN", args.run,"."

###### --- Online analysis
#ds = psana.DataSource('shmem=AMO.0:stop=no')
#for run in ds.runs():
#    for evt in run.events():
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

        # do x projection
        opal_thresh_xproj = np.sum(opal_thresh,axis=1)

        # sum up the projected image in bin range 100 to 200
        #integrated_projx = np.sum(opal_thresh_xproj[100:200])

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
        ##### Filling Histograms
        # Hit Histogram
        hithist.data = np.zeros(hithist.nbinx)
        # Hit Histogram Jitter Corrected
        hithistcorr.data = np.zeros(hithistcorr.nbinx)
        for hit in c:
            hithist.fill(hit[0])
            hithistcorr.fill(hit[0]+shift)

       
        ##### Summing over single shots (good practice to be compatible with online)
        if 'spectrum_sum' in locals():
            spectrum_sum += hithist.data
        else:
            spectrum_sum = hithist.data

        if 'spectrumcorr_sum' in locals():
            spectrumcorr_sum += hithistcorr.data
        else:
            spectrumcorr_sum = hithistcorr.data

        if 'opal_thresh_xproj_sum' in locals():
            opal_thresh_xproj_sum += opal_thresh_xproj
        else:
            opal_thresh_xproj_sum = opal_thresh_xproj

        # Debugging
        if rank==0:
            print evtGood


############### create zero arrays and dump for master
if not 'spectrum_sum_all' in locals():
    spectrum_sum_all = np.zeros_like(spectrum_sum)
comm.Reduce(spectrum_sum,spectrum_sum_all)

if not 'spectrumcorr_sum_all' in locals():
    spectrumcorr_sum_all = np.zeros_like(spectrumcorr_sum)
comm.Reduce(spectrumcorr_sum,spectrumcorr_sum_all)

if not 'opal_thresh_xproj_sum_all' in locals():
    opal_thresh_xproj_sum_all = np.zeros_like(opal_thresh_xproj_sum)
comm.Reduce(opal_thresh_xproj_sum,opal_thresh_xproj_sum_all)


evtGood_array = np.array([evtGood])
if not 'evtGood_sum_all' in locals():
    evtGood_sum_all = np.zeros_like(evtGood_array)
comm.Reduce(evtGood_array,evtGood_sum_all)


if rank==0:
    ###### Average spectra
    spectrum_sum_all_avg = spectrum_sum_all/evtGood_sum_all[0]
    spectrumcorr_sum_all_avg = spectrumcorr_sum_all/evtGood_sum_all[0]
    opal_thresh_xproj_sum_all_avg = opal_thresh_xproj_sum_all/evtGood_sum_all[0]

    ###### calculating moments of the hithist.data
    m,s = moments(opal_thresh_xproj_sum_all_avg)
    bm,bs = moments(spectrum_sum_all_avg)
    bmc,bsc = moments(spectrumcorr_sum_all_avg)

    ###### Output (debugging and info)
    print evtGood_sum_all[0], 'total events processed.' # Real number of good events
    #print evtGood*size, 'total events processed.' # Roughly ok.
    print 'Standard deviation of the raw spectrum is', round(s,1), 'pixels, the binned spectras is', round(bs,1), 'bins and the corrected, binned spectra has a standard deviation of', round(bsc,1), "bins."

    ##### Saving spectra to textfiles #####
    #names
    spectra_name = 'data/run'+str(args.run).zfill(4)+'_PixelSpectra'
    spectrabin_name = 'data/run'+str(args.run).zfill(4)+'_BinnedSpectra'
    spectrabincorr_name = 'data/run'+str(args.run).zfill(4)+'_BinnedCorrectedSpectra'
    standarddeviation_name = 'data/run'+str(args.run).zfill(4)+'_StandardDeviation_Pixel_Binned_BinnedCorrected'
    #saveing files
    np.savetxt(spectra_name,opal_thresh_xproj_sum_all_avg)
    np.savetxt(spectrabin_name,spectrum_sum_all_avg)
    np.savetxt(spectrabincorr_name,spectrumcorr_sum_all_avg)
    np.savetxt(standarddeviation_name,np.array([s,bs,bsc]))

    print 'Files have been written into the working_directory/data/.'
    
    ############### Exquisit plotting##################
    #ax = range(0,len(spectrum_sum_all_avg))
    #spectra_plot = XYPlot(evtGood_sum_all[0], "Full Run Average Spectrum", ax, spectrum_sum_all_avg,formats=".-") # make a 1D plot
    #publish.send("FULLSPECTRA", spectra_plot) # send to the display

                