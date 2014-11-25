import pyspeckit
import numpy as np
from pyspeckit.spectrum.models import inherited_voigtfitter
import matplotlib.pyplot



c = np.loadtxt('/Users/Max/ownCloud/Doktor/amoi0214_data/run0080_PixelSpectra')


#xarr = pyspeckit.spectrum.units.SpectroscopicAxis(np.linspace(-100,100,500),unit='km/s',refX=1e9,refX_units='Hz')
xarr = pyspeckit.spectrum.units.SpectroscopicAxis(range(0,len(c)),unit='Hz',refX=None,refX_units=None)


sp1 = pyspeckit.Spectrum(xarr=xarr, data=c, error=np.ones(xarr.shape[0])/20.)

sp1.plotter()

sp1.specfit(fittype='gaussian', guesses=[85000,435,10,15000,575,8,3000,950,3],
           composite_fit_color='b', clear=False, annotate=False,multifit=False)
#sp1.specfit(fittype='lorentzian', guesses=[85000,435,10],
#           composite_fit_color='g', clear=False, annotate=False)

#sp1.specfit(fittype='voigt', guesses=[85000,435,10,4,15000,575,8,3,3000,950,3,1],
#           composite_fit_color='r',clear=False,annotate=True)
#sp1.baseline(excludefit=True)
#sp1.baseline.annotate()

# this approach doesn't work right now, but it will (there's a bug I'm working on)
# it's a lot more verbose, so it's kinda ugly, but it is (in principle) more flexible
# parinfo = pyspeckit.parinfo.ParinfoList()
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='AMP',value=85000))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='FREQ',value=435))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='GWIDTH',value=10))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='LWIDTH',value=4))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='AMP',value=15000))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='FREQ',value=575))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='GWIDTH',value=8))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='LWIDTH',value=3))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='AMP',value=3000))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='FREQ',value=950))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='GWIDTH',value=3))
# parinfo.append(pyspeckit.parinfo.Parinfo(parname='LWIDTH',value=1))

# sp1.specfit(fittype='voigt', parinfo=parinfo, multifit=True,
#             composite_fit_color='c', clear=False, annotate=True)

# fwhm = sp1.specfit.fitter.analytic_fwhm()
# real_fwhm = [inherited_voigtfitter.voigt_fwhm(6.5,0.5),
#              inherited_voigtfitter.voigt_fwhm(1.5,6.5)]
# print(fwhm,real_fwhm)

matplotlib.pyplot.show()