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

def XXX(plt):
  GS = GS_ROI(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
  GS.plot_them(plt,f1="b.",f2="b.")
  GS.plot_them(plt,f1="b-",f2="b-")
  GS = GS_ROI(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  GS.plot_them(plt,f1="r.",f2="r.")
  GS.plot_them(plt,f1="r-",f2="r-")
  GS = GS_ROI(full_path("data_sherrell/Fe_fake.dat")) # with interpolated points
  GS.plot_them(plt,f1="m-",f2="m-")
  #plt.axes().set_xlim((7088,7152))
  plt.axes().set_ylim((-8.3,4.2))
  plt.show()

def restrain_II(plt):
  OX = GS_ROI(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
  RD = GS_ROI(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  MT = GS_ROI(full_path("data_sherrell/Fe_fake.dat")) # with interpolated points
  r_mean = flex.double(200)
  r_sigma = flex.double(200)
  for ichannel in range(100):
    fp_pop = flex.mean_and_variance(flex.double([OX.fp[ichannel],RD.fp[ichannel],MT.fp[ichannel]]))
    fdp_pop = flex.mean_and_variance(flex.double([OX.fdp[ichannel],RD.fdp[ichannel],MT.fdp[ichannel]]))
    r_mean[ichannel] = fp_pop.mean(); r_mean[100+ichannel] = fdp_pop.mean()
    r_sigma[ichannel] = fp_pop.unweighted_sample_standard_deviation()
    r_sigma[100+ichannel] = fdp_pop.unweighted_sample_standard_deviation()

  plt.stackplot(OX.energy,r_mean[0:100]-r_sigma[0:100],color=('lightgreen'))
  plt.stackplot(OX.energy,r_mean[0:100]+r_sigma[0:100],color=('white'))
  plt.stackplot(OX.energy,r_mean[100:200]+r_sigma[100:200],color=('lightgreen'))
  plt.stackplot(OX.energy,r_mean[100:200]-r_sigma[100:200],color=('white'))
  plt.plot(OX.energy,r_mean[0:100]+r_sigma[0:100],'g-')
  plt.plot(OX.energy,r_mean[100:200]+r_sigma[100:200],'g-')
  plt.plot(OX.energy,r_mean[0:100]-r_sigma[0:100],'g-')
  plt.plot(OX.energy,r_mean[100:200]-r_sigma[100:200],'g-')
  plt.plot(OX.energy,r_mean[0:100],'g-')
  plt.plot(OX.energy,r_mean[100:200],'g-')
  OX.plot_them(plt,f1="b.",f2="b.")
  OX.plot_them(plt,f1="b-",f2="b-")
  RD = GS_ROI(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  RD.plot_them(plt,f1="r.",f2="r.")
  RD.plot_them(plt,f1="r-",f2="r-")
  MT = GS_ROI(full_path("data_sherrell/Fe_fake.dat")) # with interpolated points
  MT.plot_them(plt,f1="m-",f2="m-")
  plt.axes().set_ylim((-8.3,4.2))
  plt.show()

def restrain_II_values():
  OX = GS_ROI(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
  RD = GS_ROI(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  MT = GS_ROI(full_path("data_sherrell/Fe_fake.dat")) # with interpolated points
  r_mean = flex.double(200)
  r_sigma = flex.double(200)
  for ichannel in range(100):
    fp_pop = flex.mean_and_variance(flex.double([OX.fp[ichannel],RD.fp[ichannel],MT.fp[ichannel]]))
    fdp_pop = flex.mean_and_variance(flex.double([OX.fdp[ichannel],RD.fdp[ichannel],MT.fdp[ichannel]]))
    r_mean[ichannel] = fp_pop.mean(); r_mean[100+ichannel] = fdp_pop.mean()
    r_sigma[ichannel] = fp_pop.unweighted_sample_standard_deviation()
    r_sigma[100+ichannel] = fdp_pop.unweighted_sample_standard_deviation()
  return r_mean, r_sigma

restrain_II_r_mean, restrain_II_r_sigma = restrain_II_values()

def restrain_II_compute_functional_and_gradients(values, mean, sigma):
  # for example: values[400] for fp,fdp of FE1,FE2; mean[200] for fp,fdp; sigma[200] for fp,fdp
  f = 0.
  g1 = flex.double(len(mean))
  g2 = flex.double(len(mean))
  for idx in range(200):
    f += (0.5/(sigma[idx]*sigma[idx]) * (values[idx] - mean[idx])**2)
    f += (0.5/(sigma[idx]*sigma[idx]) * (values[200+idx] - mean[idx])**2)
  g1 = (values[0:200]-mean)/(sigma*sigma)
  g2 = (values[200:400]-mean)/(sigma*sigma)
  return f,g1,g2

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
  print("mean",mean,"sigma",sigma)
  compute_functional_and_gradients_fp(FE1_fp=OX.fp,FE2_fp=RD.fp,mean=mean,sigma=sigma)
  n,bins,patches = plt.hist(delta_fp, 50, normed=0, facecolor="orange", alpha=0.75)
  plt.plot(plotx, ploty, "r-")
  plt.xlabel("Delta fp")
  plt.title("Histogram of Delta fp")
  plt.axis([-1,1,0,40])
  plt.show()

def compute_functional_and_gradients_fp(FE1_fp,FE2_fp,mean,sigma,constrain_endpoints=False):
  f = 0.
  if constrain_endpoints:  effective_range=range(1,99)
  else: effective_range=range(100)
  for i in effective_range:
    f += (0.5/(sigma*sigma) * (FE1_fp[i] - FE1_fp[i-1] -mean)**2)
    f += (0.5/(sigma*sigma) * (FE2_fp[i] - FE2_fp[i-1] -mean)**2)
  g1 = flex.double(len(FE1_fp))
  g2 = flex.double(len(FE2_fp))
  for i in range(100):
    g1[i] += 2. * FE1_fp[i]; g2[i] += 2. * FE2_fp[i]
    if i > 0:
      g1[i] -= FE1_fp[i-1]; g2[i] -= FE2_fp[i-1]
    else:
      g1[i] -= FE1_fp[i]; g2[i] -= FE2_fp[i]
    if i < 99:
      g1[i] -= FE1_fp[i+1]; g2[i] -= FE2_fp[i+1]
    else:
      g1[i] -= FE1_fp[i]; g2[i] -= FE2_fp[i]

  g1 /= sigma*sigma; g2/= sigma*sigma
  return f,g1,g2

def tst_analytical_fp():
  from libtbx.test_utils import approx_equal
  OX = GS_ROI(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
  RD = GS_ROI(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  fbase,g1base,g2base = compute_functional_and_gradients_fp(
         FE1_fp=OX.fp,FE2_fp=RD.fp,mean=0.0,sigma=0.1)
  EPS = 0.0001
  for ichannel in range(100):
    OX.fp[ichannel]+=EPS
    fnew,g1new,g2new = compute_functional_and_gradients_fp(
         FE1_fp=OX.fp,FE2_fp=RD.fp,mean=0.0,sigma=0.1)
    assert approx_equal((fnew-fbase)/EPS , g1base[ichannel],eps = 1e-1)
    OX.fp[ichannel]-=EPS
    RD.fp[ichannel]+=EPS
    fnew,g1new,g2new = compute_functional_and_gradients_fp(
         FE1_fp=OX.fp,FE2_fp=RD.fp,mean=0.0,sigma=0.1)
    assert approx_equal((fnew-fbase)/EPS , g2base[ichannel],eps = 1e-1)
    RD.fp[ichannel]-=EPS

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
  ploty = flex.double([displaynorm *
       (1./math.sqrt(2.*math.pi*sigma*sigma))*math.exp(-0.5*(math.pow(x-mean,2))/(sigma*sigma))
       for x in plotx])
  print("mean",mean,"sigma",sigma)
  #compute_functional_and_gradients_fp(FE1_fp=OX.fp,FE2_fp=RD.fp,mean=mean,sigma=sigma)
  n,bins,patches = plt.hist(delta_fdp, 50, normed=0, facecolor="orange", alpha=0.75)
  plt.plot(plotx, ploty, "r-")
  plt.xlabel("Delta d")
  plt.title("Histogram of Delta fdp")
  plt.axis([-1,1,0,40])
  plt.show()

if __name__=="__main__":

  from matplotlib import pyplot as plt
  #XXX(plt)
  restrain_II(plt)
  #fp_distro(plt)
  tst_analytical_fp()
  #fdp_distro(plt)
