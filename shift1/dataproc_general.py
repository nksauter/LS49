from __future__ import division
from six.moves import range
#source /reg/g/psdm/etc/psconda.sh
from psana import *
import numpy as np
import sys
import math
proposal = sys.argv[1]

def gen_ls49(ds):
  for run in ds.runs():
    print run
    times = run.times()
    for nevent, t in enumerate(times):
      evt = run.event(t)
      yield nevent,evt,t
def gen_lg36(ds):
  for nevent,evt in enumerate(ds.events()):
    yield nevent,evt,None

from LS49.spectra.energy_results import linear_fit
class calib_A(linear_fit):
  def detidx_as_ebeam(self,x):
    return self.m*x+self.c
  def ebeam_as_detidx(self,y):
    return (1./self.m)*y - (self.c/self.m)
  def plot(self,plt):
    print "expidx vs ebeam energy"
    e1 = self.detidx_as_ebeam(200.)
    e2 = self.detidx_as_ebeam(1800.)
    idx1 = self.ebeam_as_detidx(e1)
    idx2 = self.ebeam_as_detidx(e2)
    plt.plot(self.y, self.x, "b.")
    plt.plot([e1,e2],[idx1,idx2],"r-")
    plt.show()

class GaussFit:
  def __init__(self,data):
    self.data = data
    numerator = np.sum( range(len(self.data)) * self.data )
    denominator = np.sum( self.data )
    self.mean_index = numerator/denominator
    sqnum = np.zeros(len(self.data))
    for ix in range(len(self.data)):
      sqnum[ix] = (float(ix) - self.mean_index)**2 * self.data[ix]
    numerator = np.sum( sqnum )
    self.stddevidx = math.sqrt(numerator/denominator)



sources = dict(lg36={"src":'exp=cxig3614:run=209',"fee":'Fee_Orca_Spectrometer',"gen":gen_lg36},
               ls49={"src":"exp=mfxls4916:run=23:idx","fee":"XrayTransportDiagnostic.0:OrcaFl40.0","gen":gen_ls49}
              )

ds = DataSource(sources[proposal]["src"]) #'exp=cxig3614:run=209')
#det = Detector('Ds1CsPad')
det = Detector(sources[proposal]["fee"]) #'Fee_Orca_Spectrometer')
ebeamDet = Detector('EBeam')
energy = []
expidx = []
spectra = []
with_energy=0
plots = 3

for nevent,evt,time in sources[proposal]["gen"](ds):
    if with_energy==1010:break
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
    print "Nevent %d, energy %8.2feV, w/energy=%d"%(nevent, ebeam.ebeamPhotonEnergy(),with_energy)

    import matplotlib.pyplot as plt
    #plt.imshow(img,vmin=-2,vmax=2) # normalization of luminance (doesnt work)
    #plt.show()
    if plots>5:
      plt.imshow(img, interpolation='nearest')
      plt.show()

    summed  = img.sum(axis=0)
    lower = summed[0:50].mean()
    upper = summed[-50:].mean()
    baseline = (np.array(range(len(summed)))/len(summed))*(upper-lower)+lower # baseline uses the left & right edges
    min_summed = np.min(summed)
    #print "minimum ",min_summed
    real = np.array(list(summed-baseline))
    if plots>4:
      plt.plot(range(len(real)), real, 'r-')
      plt.show() # show baseline-corrected spectrum

    fr = np.fft.rfft(real)
    #print type(fr), len(real)//2, len(fr.real)
    if plots>5:
      plt.plot(range(len(fr.real)), fr.real, 'b-')
      plt.show()
    #low_pass_fr: # highest 3/4 are zeroed out
    for x in range(len(fr)//4, len(fr)):
      fr[x]=0.+0.j
    filtered_real = np.fft.irfft(fr)
    # get the expectation value of the index
    numerator = np.sum( range(len(filtered_real)) * filtered_real )
    denominator = np.sum( filtered_real )
    mean_index = numerator/denominator
    maxplot = np.max(filtered_real)
    if plots>3:
      #plt.plot(xrange(len(real)), real, 'r-')
      plt.plot(range(len(filtered_real)), filtered_real, 'g-')
      plt.plot([mean_index,mean_index],[0.15*maxplot, 1.05*maxplot],'k-')
      plt.show()
    expidx.append(mean_index)
    energy.append(ebeam.ebeamPhotonEnergy())
    spectra.append(filtered_real.astype('float32'))

# plotting the correlation between ebeam (x) and mean spectrometer index (y)
CA = calib_A(data = dict(expidx=expidx,energy=energy))
#CA.plot(plt)

# plot average spectrum
nchannel = len(spectra[0])
sumspectra = np.zeros(nchannel)
for idsp in range(len(spectra)):
  sumspectra+=spectra[idsp]
#plt.plot(xrange(nchannel),sumspectra,"r-") # raw spectrometer channel, no calib
xenergy=np.zeros(nchannel)
for ix in range(nchannel):
  xenergy[ix] = CA.detidx_as_ebeam(ix)
statsA = GaussFit(sumspectra)
#print statsA.mean_index, statsA.stddevidx
meanE = CA.detidx_as_ebeam(statsA.mean_index)
fullheight = sumspectra[int(statsA.mean_index)]
plt.plot([meanE,meanE],[0.,fullheight],'r-')
print "Mean energy %8.2feV"%meanE
lower = CA.detidx_as_ebeam(statsA.mean_index- statsA.stddevidx)
upper = CA.detidx_as_ebeam(statsA.mean_index+ statsA.stddevidx)
print "Full width %8.2feV"%(upper-lower)
plt.plot([lower,upper],[fullheight/2.,fullheight/2.],'r-')

plt.plot(xenergy,sumspectra,"b-")
plt.show()
