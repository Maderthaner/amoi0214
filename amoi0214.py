
# useful web pages:
# psana-python: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual
# realtime plots: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-Real-timeOnlinePlotting/Monitoring
# mpi parallelization: https://confluence.slac.stanford.edu/display/PSDM/psana+-+Python+Script+Analysis+Manual#psana-PythonScriptAnalysisManual-MPIParallelization

# plotting commands (when run on the same node as the analysis):
# psplot OPAL &
# psplot XPROJ &

# plotting commands (when run on the a different node) '-s' means 'server':
# psplot -s daq-amo-mon03 OPAL &




# Standard PYTHON modules
print "IMPORTING STANDARD PYTHON MODULES...",
import matplotlib.pyplot as plt
import numpy as np
import math
import collections
print "DONE"

# LCLS psana to read data
print "IMPORTING psana...",
import psana 
print "DONE"

# For online plotting
from psmon import publish
from psmon.plots import XYPlot,Image

# custom algorithms
from pypsalg import find_blobs

def moments(arr):
    if np.count_nonzero(arr) == 0:
        return 0,0
    nbins = len(arr)
    bins  = range(nbins)
    mean  = np.average(bins,weights=arr)
    var   = np.average((bins-mean)**2,weights=arr)
    sdev  = np.sqrt(var)
    return mean,sdev

publish.init()

# --- Offline analysis
ds = psana.DataSource("exp=AMO/amoi0214:run=1")    

# --- Online analysis
# = psana.DataSource('shmem=AMO.0:stop=no')

# Set up detectors to read
#iTOF_Src = psana.Source("DetInfo(AmoITOF.0:Acqiris.0)")
opal_src = psana.Source("DetInfo(AmoEndstation.1:Opal1000.0)")

# Set up event counters
eventCounter = 0
evtGood = 0
evtBad = 0
evtUsed = 0

# Buffers for histories
opal_hit_buff = collections.deque(maxlen=25)
opal_hit_avg_buff = collections.deque(maxlen=25)
opal_circ_buff = collections.deque(maxlen=25)
xproj_int_buff = collections.deque(maxlen=25)
xproj_int_avg_buff = collections.deque(maxlen=25)
moments_buff = collections.deque(maxlen=25)
moments_avg_buff = collections.deque(maxlen=25)



# Loop for online/offline (no counter)
for evt in ds.events():
# Loop for online/offline (with counter)
#for eventCounter, evt in enumerate(ds.events()) :

    opal_raw = evt.get(psana.Camera.FrameV1,opal_src)

    # Check all detectors are read in
    if (opal_raw is None) :
        evtBad += 1
        continue
    else:
        evtGood += 1

    ####### Copy detector data for further analysis #########
    opal = np.copy(opal_raw.data16())

    # threshold the image
    threshold = 1
    opal_thresh = (opal>threshold)*opal

    # do two projections
    opal_thresh_xproj = np.sum(opal_thresh,axis=1)
    opal_thresh_yproj = np.sum(opal_thresh,axis=0)

    # calculating moments of the xproj
    m,s = moments(opal_thresh_xproj)


    # accumulate sums
    # if 'xproj_sum' in locals():
    #     xproj_sum += opal_thresh_xproj
    # else:
    #     xproj_sum = opal_thresh_xproj
    # if 'thresh_sum' in locals():
    #     thresh_sum += opal_thresh
    # else:
    #     thresh_sum = opal_thresh

    # # compute average x-projection
    # xproj_avg = xproj_sum/evtGood

    #plt.plot(opal_thresh_xproj)
    #plt.show()

    print evtGood

    # sum up the projected image in bin range 100 to 200
    integrated_projx = np.sum(opal_thresh_xproj[100:200])
    

    # do blob finding
    thresh_sigma = 5
    c,w = find_blobs.find_blobs(opal,thresh_sigma)
    print 'number of hits:',len(c)

    ############## Histories of certain values #############
    # Opal history
    opal_circ_buff.append(opal)
    opal_sum = sum(opal_circ_buff)

    # Opal hitcounter history
    opal_hit_buff.append(len(c))
    opal_hit_avg = sum(opal_hit_buff)/len(opal_hit_buff)
    opal_hit_avg_buff.append(opal_hit_avg)

    # x-projection history
    xproj_int_buff.append(integrated_projx)
    xproj_int_avg = sum(xproj_int_buff)/len(xproj_int_buff)
    xproj_int_avg_buff.append(xproj_int_avg)

    # moments history
    moments_buff.append(s)
    moments_avg = sum(moments_buff)/len(moments_buff)
    moments_avg_buff.append(moments_avg)
    


    ############# Things that are done after the hitfinder ############
    if len(c)==0:
        print '****** no hits ******'
        continue

    #print 'moments:',moments(xproj_avg)
    

    #find_blobs.draw_blobs(opal,c,w) # draw the blobs in the opal picture

    # only update the plots "once in a while"
    if evtGood%10 == 0:
        ax = range(0,len(opal_thresh_xproj))
        xproj_plot = XYPlot(evtGood, "X Projection", ax, opal_thresh_xproj) # make a 1D plot
        publish.send("XPROJ", xproj_plot) # send to the display
        #
        ax = range(0,len(opal_thresh_yproj))
        yproj_plot = XYPlot(evtGood, "Y Projection", ax, opal_thresh_yproj) # make a 1D plot
        publish.send("YPROJ", yproj_plot) # send to the display
        #
        opal_hit_avg_arr = np.array(list(opal_hit_avg_buff))
        ax2 = range(0,len(opal_hit_avg_arr))
        hitrate_avg_plot = XYPlot(evtGood, "Hitrate history", ax2, opal_hit_avg_arr) # make a 1D plot
        publish.send("HITRATE", hitrate_avg_plot) # send to the display
        #
        xproj_int_avg_arr = np.array(list(xproj_int_avg_buff))
        ax3 = range(0,len(xproj_int_avg_arr))
        xproj_int_plot = XYPlot(evtGood, "XProjection running avg history", ax3, xproj_int_avg_arr) # make a 1D plot
        publish.send("XPROJRUNAVG", xproj_int_plot) # send to the display
        #
        moments_avg_arr = np.array(list(moments_avg_buff))
        ax4 = range(0,len(moments_avg_arr))
        moments_avg_plot = XYPlot(evtGood, "Moments running avg history", ax4, moments_avg_arr) # make a 1D plot
        publish.send("MOMENTSRUNAVG", moments_avg_plot) # send to the display
        #
        opal_plot = Image(evtGood, "Opal", opal) # make a 2D plot
        publish.send("OPAL", opal_plot) # send to the display
        #
        opal_plot_avg = Image(evtGood, "Opal Average", opal_sum) # make a 2D plot
        publish.send("OPALAVG", opal_plot_avg) # send to the display

        