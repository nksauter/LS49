from __future__ import division, absolute_import

from LS49.sim.fdp_plot import george_sherrell

if __name__=="__main__":

  from matplotlib import pyplot as plt

  GS = george_sherrell("data_sherrell/pf-rd-ox_fftkk.out")
  GS.plot_them(plt,f1="b.",f2="b.")
  GS.plot_them(plt,f1="b-",f2="b-")
  GS = george_sherrell("data_sherrell/pf-rd-red_fftkk.out")
  GS.plot_them(plt,f1="r.",f2="r.")
  GS.plot_them(plt,f1="r-",f2="r-")
  GS = george_sherrell("data_sherrell/Fe.dat")
  GS.plot_them(plt,f1="m-",f2="m-")
  plt.axes().set_xlim((7088,7152))
  plt.axes().set_ylim((-8.3,4.2))


  print list(GS.energy)
  print list(GS.fdp)

  plt.show()
