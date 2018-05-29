from __future__ import division
#source /reg/g/psdm/etc/psconda.sh
from psana import *
import numpy as np
ds = DataSource('exp=cxig3614:run=209')
#det = Detector('Ds1CsPad')
det = Detector('Fee_Orca_Spectrometer')
ebeamDet = Detector('EBeam')
energy = []
expidx = []
spectra = []
with_energy=0
for nevent,evt in enumerate(ds.events()):
    if with_energy==101010:break
    # includes pedestal subtraction, common-mode correction, bad-pixel
    # suppresion, and returns an "unassembled" 3D array of cspad panels
    calib_array = det.calib(evt)
    # this is the same as the above, but also uses geometry to
    # create an "assembled" 2D image (including "fake pixels" in gaps)
    img = det.image(evt)
    if img is None:
      print 'None',nevent
      continue
    ebeam = ebeamDet.get(evt)
    if ebeam is None:
      print 'None ebeam',nevent
      continue
    with_energy+=1
    print nevent, ebeam.ebeamPhotonEnergy(),with_energy

    import matplotlib.pyplot as plt
    #from IPython import embed; embed()
    #plt.imshow(img,vmin=-2,vmax=2)
    ##plt.imshow(img, interpolation='nearest')
    ##plt.show()
    #from IPython import embed; embed()
    summed  = img.sum(axis=0)
    lower = summed[0:50].mean()
    upper = summed[-50:].mean()
    baseline = (np.array(xrange(len(summed)))/len(summed))*(upper-lower)+lower
    min_summed = np.min(summed)
    print "minimum ",min_summed
    real = np.array(list(summed-baseline))
    #from IPython import embed; embed()
    #plt.plot(xrange(len(real)), real, 'r-')
    #plt.show()
    fr = np.fft.rfft(real)
    print type(fr), len(real)/2, len(fr.real)
    #plt.plot(xrange(len(fr.real)), fr.real, 'b-')
    #plt.show()
    #low_pass_fr:
    for x in xrange(len(fr)/4, len(fr)):
      fr[x]=0.+0.j
    filtered_real = np.fft.irfft(fr)
    # get the expectation value of the index
    numerator = np.sum( xrange(len(filtered_real)) * filtered_real )
    denominator = np.sum( filtered_real )
    mean_index = numerator/denominator
    maxplot = np.max(filtered_real)
    #plt.plot(xrange(len(real)), real, 'r-')
    #plt.plot(xrange(len(filtered_real)), filtered_real, 'g-')
    #plt.plot([mean_index,mean_index],[0.15*maxplot, 1.05*maxplot],'k-')
    expidx.append(mean_index)
    energy.append(ebeam.ebeamPhotonEnergy())
    spectra.append(filtered_real.astype('float32'))
    #plt.show()
import pickle
pickle.dump(dict(expidx=expidx,energy=energy,spectra=spectra), open("/reg/d/psdm/cxi/cxig3614/scratch/spectra209_101010.pickle","wb"), pickle.HIGHEST_PROTOCOL)
#plt.plot(energy, expidx, "b.")
#plt.show()
