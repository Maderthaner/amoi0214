import psana 
import numpy as np
ds = psana.DataSource("exp=AMO/amoi0214:run=2")
opal_src = psana.Source("DetInfo(AmoEndstation.1:Opal1000.0)")
evtBad = 0
evtGood = 0
for evt in ds.events():
    opal_raw = evt.get(psana.Camera.FrameV1,opal_src)
    if (opal_raw is None) :
        evtBad += 1
        continue
    else:
        evtGood += 1
    if not 'opalsum' in locals():
        opalsum = opal_raw.data16().astype(float)
    else:
        opalsum += opal_raw.data16()
    if evtGood>99:
        break
opalped = opalsum/float(evtGood)
np.save('ped.npy',opalped)
print 'Wrote pedestals to ped.npy, averaging %d events.' %evtGood
