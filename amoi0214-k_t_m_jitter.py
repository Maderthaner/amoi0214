
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
opal_src = psana.Source("DetInfo(AmoEndstation.0:Opal1000.0)")
acq_src = psana.Source("DetInfo(AmoITOF.0:Acqiris.0)")

# Set up event counters
eventCounter = 0
evtGood = 0
evtBad = 0
evtUsed = 0
evtSeeded = 0

# Buffers for histories
history_len = 6
history_len_long = 50
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


from Histogram import hist1d

hithist = hist1d(100,0.,1024.)
hitjitter = hist1d(100,0.,1024.)

# ion yield array for sxrss scan
ion_yield = np.zeros(102) ##choose appropriate range

###### --- Online analysis
ds = psana.DataSource('shmem=AMO.0:stop=no')
for run in ds.runs():
    for evt in run.events():
###### --- End of Online analysis region

###### --- Offline analysis
#ds = psana.DataSource("exp=AMO/amoi0214:run=2:idx")
#for run in ds.runs():
#    times = run.times()
#    mylength = len(times)/size
#    mytimes= times[rank*mylength:(rank+1)*mylength]
#    for i in range(mylength):
#        evt = run.event(mytimes[i])
###### --- End of Offline analysis region

        opal_raw = evt.get(psana.Camera.FrameV1,opal_src)

        acq = evt.get(psana.Acqiris.DataDescV1,acq_src)
        acqnumchannels = acq.data_shape()
        chan=0

        wf = acq.data(chan).waveforms()[0]
        
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
            #print blob_centroid
            shift = 512-int(blob_centroid)
            index = int(blob_centroid/100)        
            mu,sig = moments(x2)	
	
        else:
            shift = 0
            
        #print shift
        # Hit Histogram
        hithist.data = np.zeros(hithist.nbinx)
        hitjitter.data = np.zeros(hitjitter.nbinx)
        for hit in c:
            hithist.fill(hit[0]+shift)
            hitjitter.fill(hit[0])
            

        #print c

        #print 'number of hits:',len(c)
        #find_blobs.draw_blobs(opal,c,w) # draw the blobs in the opal picture

        # find yield of waveform over certain interval
        wf_yield = np.sum(wf[0:60]) #choose proper window later

        if 'sig' in locals():
            if sig < 10: #change width condition
                evtSeeded += 1
                ion_yield[index] = ion_yield[index] + wf_yield
   

        ############## Histories of certain values #############

        # x-projection histogram history
        hitxprojhist_buff.append(hithist.data)
        hitxprojhist_sum = sum(hitxprojhist_buff)
        hitxprojjitter_buff.append(hitjitter.data)
        hitxprojjitter_sum = sum(hitxprojjitter_buff)


        # Opal hitcounter history
        opal_hit_buff.append(len(c))
        opal_hit_sum = np.array([sum(opal_hit_buff)])#/len(opal_hit_buff)
        

        # x-projection history
        #xproj_int_buff.append(integrated_projx)
        #xproj_int_sum = np.array([sum(xproj_int_buff)])#/len(xproj_int_buff)
        

        # Opal history
        opal_circ_buff.append(opal)
        opal_sum = sum(opal_circ_buff)
        

        # only update the plots and call conm.Reduce "once in a while"
        if evtGood%6 == 0:
            ### create empty arrays and dump for master
            #if not 'moments_sum_all' in locals():
            #    moments_sum_all = np.empty_like(moments_sum)
            #comm.Reduce(moments_sum,moments_sum_all)

            #if not 'xproj_in_sum_all' in locals():
            #    xproj_int_sum_all = np.empty_like(xproj_int_sum)
            #comm.Reduce(xproj_int_sum,xproj_int_sum_all)

            if not 'opal_hit_sum_all' in locals():
                opal_hit_sum_all = np.empty_like(opal_hit_sum)
            comm.Reduce(opal_hit_sum,opal_hit_sum_all)

            if not 'opal_sum_all' in locals():
                opal_sum_all = np.empty_like(opal_sum)
            comm.Reduce(opal_sum,opal_sum_all)

            if not 'hithist_sum_all' in locals():
                hithist_sum_all = np.empty_like(hithist.data)
            comm.Reduce(hithist.data,hithist_sum_all)

            if not 'hitxprojhist_sum_all' in locals():
                hitxprojhist_sum_all = np.empty_like(hitxprojhist_sum)
            comm.Reduce(hitxprojhist_sum,hitxprojhist_sum_all)

            if not 'hitjitter_sum_all' in locals():
                hitjitter_sum_all = np.empty_like(hitjitter.data)
            comm.Reduce(hitjitter.data,hitjitter_sum_all)

            if not 'hitxprojjitter_sum_all' in locals():
                hitxprojjitter_sum_all = np.empty_like(hitxprojjitter_sum)
            comm.Reduce(hitxprojjitter_sum,hitxprojjitter_sum_all)

            if not 'ion_yield_sum_all' in locals():
                ion_yield_sum_all = np.empty_like(ion_yield)
            comm.Reduce(ion_yield,ion_yield_sum_all)


            if rank==0:
                ###### calculating moments of the hithist.data
                m,s = moments(hitxprojjitter_sum_all)

                ###### History on master
                print eventCounter*size, 'total events processed.'
                opal_hit_avg_buff.append(opal_hit_sum_all[0]/(len(opal_hit_buff)*size))
                print evtSeeded/evtGood 'is the ration of seeded/all shots'
                #xproj_int_avg_buff.append(xproj_int_sum_all[0]/(len(xproj_int_buff)*size))
                #moments_avg_buff.append(moments_sum_all[0]/(len(moments_buff)*size))
                ### moments history
                moments_buff.append(s)
                #moments_sum = np.array([sum(moments_buff)])#/len(moments_buff)

                # Exquisit plotting
                ax = range(0,len(opal_thresh_xproj))
                xproj_plot = XYPlot(evtGood, "X Projection", ax, opal_thresh_xproj,formats=".-") # make a 1D plot
                publish.send("XPROJ", xproj_plot) # send to the display

                # #
                opal_hit_avg_arr = np.array(list(opal_hit_avg_buff))
                # print opal_hit_avg_buff
                ax2 = range(0,len(opal_hit_avg_arr))
                hitrate_avg_plot = XYPlot(evtGood, "Hitrate history", ax2, opal_hit_avg_arr,formats=".-") # make a 1D plot
                # print 'publish',opal_hit_avg_arr
                publish.send("HITRATE", hitrate_avg_plot) # send to the display
                # #
                #xproj_int_avg_arr = np.array(list(xproj_int_avg_buff))
                #ax3 = range(0,len(xproj_int_avg_arr))
                #xproj_int_plot = XYPlot(evtGood, "XProjection running avg history", ax3, xproj_int_avg_arr) # make a 1D plot
                #publish.send("XPROJRUNAVG", xproj_int_plot) # send to the display
                # #
                moments_avg_arr = np.array(list(moments_buff))
                ax4 = range(0,len(moments_avg_arr))
                moments_avg_plot = XYPlot(evtGood, "Standard Deviation of avg Histogram", ax4, moments_avg_arr,formats=".-") # make a 1D plot
                publish.send("MOMENTSRUNAVG", moments_avg_plot) # send to the display

                # average histogram of hits
                ax5 = range(0,hitxprojhist_sum_all.shape[0])
                hithist_plot = XYPlot(evtGood, "Hit Hist2", ax5, hitxprojhist_sum_all/(len(hitxprojhist_buff)*size),formats=".-") # make a 1D plot
                publish.send("HITHIST", hithist_plot) # send to the display

                # histogram of single shot
                ax6 = range(0,hithist.data.shape[0])
                hit1shot_plot = XYPlot(evtGood, "Hit Hist 1 shot", ax6, hithist.data,formats=".-") # make a 1D plot
                publish.send("HIT1SHOT", hit1shot_plot) # send to the display

                # #
                opal_plot = Image(evtGood, "Opal", opal) # make a 2D plot
                publish.send("OPAL", opal_plot) # send to the display

                # #
                opal_plot_avg = Image(evtGood, "Opal Average", opal_sum_all/(len(opal_circ_buff)*size)) # make a 2D plot
                publish.send("OPALAVG", opal_plot_avg) # send to the display
                
                # Exquisit plotting, acqiris
                ax_wf = range(0,len(wf))
                wf_plot = XYPlot(evtGood, "TOF Trace", ax_wf, wf,formats=".-") # make a 1D plot
                publish.send("WF", wf_plot) # send to the display

                # Exquisit plotting, ion yield
                ax_ion = range(0,len(ion_yield))
                ion_yield_plot = XYPlot(evtGood, "Ion Yield", ax_ion, ion_yield,formats=".-") # make a 1D plot
                publish.send("ION", ion_yield_plot) # send to the display

                # average of ion yield sum all
                ax_ion_all = range(0,ion_yield_sum_all.shape[0])
                ion_all_plot = XYPlot(evtGood, "Ion Yield All", ax_ion_all, ion_yield_sum_all,formats=".-") # make a 1D plot
                publish.send("IONALL", ion_all_plot) # send to the display

                # histogram of single shot, no jitter correction
                axj = range(0,hitjitter.data.shape[0])
                hitjitter_plot = XYPlot(evtGood, "Hit Hist No Jitter Correction", axj, hitjitter.data,formats=".-") # make a 1D plot
                publish.send("HITJITTER", hitjitter_plot) # send to the display

               # average histogram of hits, no jitter correction
                axja = range(0,hitxprojjitter_sum_all.shape[0])
                hitjitter_plot_avg = XYPlot(evtGood, "Hit Average, No Jitter Correction", axja, hitxprojjitter_sum_all/(len(hitxprojjitter_buff)*size),formats=".-") # make a 1D plot
                publish.send("HITJITAVG", hitjitter_plot_avg) # send to the display
