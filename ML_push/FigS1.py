from __future__ import division, absolute_import, print_function
from scitbx.array_family import flex
from LS49.sim.fdp_plot import george_sherrell
import math
# %%% boilerplate context: specialize to packaged big data %%%
import os
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
from LS49.sim import step5_pad
step5_pad.big_data = ls49_big_data
from LS49.sim.step5_pad import full_path
# %%%%%%

class GS_ROI(george_sherrell):
  def __init__(self,path):
    george_sherrell.__init__(self,path)
    selection = (self.energy >= 7070.).__and__(self.energy < 7170.)
    self.energy = self.energy.select(selection)
    self.fp = self.fp.select(selection)
    self.fdp = self.fdp.select(selection)


def fp_distro(plt):
  OX = GS_ROI(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
  RD = GS_ROI(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  MM = GS_ROI(full_path("data_sherrell/Fe_fake.dat"))
  delta_fp = flex.double()

  for i in range(1,100):
    #delta_fp.append ( MM.fp[i]-MM.fp[i-1] );
    delta_fp.append ( OX.fp[i]-OX.fp[i-1] ); delta_fp.append ( RD.fp[i]-RD.fp[i-1] )
  STATS = flex.mean_and_variance(delta_fp)
  mean = STATS.mean()
  sigma = STATS.unweighted_sample_standard_deviation()
  displaynorm = 3.0
  plotx = flex.double([0.005*x for x in range(-100,100)])
  ploty = flex.double([displaynorm *
       (1./math.sqrt(2.*math.pi*sigma*sigma))*math.exp(-0.5*(math.pow(x-mean,2))/(sigma*sigma))
       for x in plotx])
  mean2 = 0.0
  sigma2 = 0.1
  displaynorm2 = 20. * math.sqrt(2 * math.pi * sigma2 *sigma2)
  ploty2 = flex.double([displaynorm2 *
       (1./math.sqrt(2.*math.pi*sigma2*sigma2))*math.exp(-0.5*(math.pow(x-mean2,2))/(sigma2*sigma2))
       for x in plotx])
  print("mean",mean,"sigma",sigma)
  #compute_functional_and_gradients_fp(FE1_fp=OX.fp,FE2_fp=RD.fp,mean=mean,sigma=sigma)
  n,bins,patches = plt.hist(delta_fp, 50, normed=0, facecolor="orange", alpha=0.75)
  plt.plot(plotx, ploty, "r-")
  plt.plot(plotx, ploty2, "b-")

  #plt.xlabel("Delta fp")
  #plt.title("Histogram of Delta fp")
  plt.axis([-0.75,0.75,0,40])


def fdp_distro(plt):
  OX = GS_ROI(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
  RD = GS_ROI(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  MM = GS_ROI(full_path("data_sherrell/Fe_fake.dat"))
  delta_fdp = flex.double()

  for i in range(1,100):
    #delta_fp.append ( MM.fp[i]-MM.fp[i-1] );
    delta_fdp.append ( OX.fdp[i]-OX.fdp[i-1] ); delta_fdp.append ( RD.fdp[i]-RD.fdp[i-1] )
  STATS = flex.mean_and_variance(delta_fdp)
  mean = STATS.mean()
  sigma = STATS.unweighted_sample_standard_deviation()
  displaynorm = 3.0
  plotx = flex.double([0.005*x for x in range(-100,100)])
  mean2 = 0.0
  sigma2 = 0.2
  displaynorm2 = 20. * math.sqrt(2 * math.pi * sigma2 *sigma2)
  ploty = flex.double([displaynorm *
       (1./math.sqrt(2.*math.pi*sigma*sigma))*math.exp(-0.5*(math.pow(x-mean,2))/(sigma*sigma))
       for x in plotx])
  ploty2 = flex.double([displaynorm2 *
       (1./math.sqrt(2.*math.pi*sigma2*sigma2))*math.exp(-0.5*(math.pow(x-mean2,2))/(sigma2*sigma2))
       for x in plotx])
  print("mean",mean,"sigma",sigma)
  #compute_functional_and_gradients_fp(FE1_fp=OX.fp,FE2_fp=RD.fp,mean=mean,sigma=sigma)
  n,bins,patches = plt.hist(delta_fdp, 50, normed=0, facecolor="orange", alpha=0.75)
  plt.plot(plotx, ploty, "r-")
  plt.plot(plotx, ploty2, "b-")
  #plt.xlabel("Delta d")
  #plt.title("Histogram of Delta fdp")
  plt.axis([-0.75,0.75,0,40])

if __name__=="__main__":

  from matplotlib import pyplot as plt
  fig, axes = plt.subplots(2, 1, figsize=(3.2,4.8)) #,figsize=(24,14)) default 6.4,4.8; new 4.8,4.8
  fp_distro(axes[0])
  fdp_distro(axes[1])
  plt.show()
