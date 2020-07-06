from __future__ import division, absolute_import, print_function

from LS49.sim.fdp_plot import george_sherrell
# %%% boilerplate context: specialize to packaged big data %%%
from LS49 import ls49_big_data
from LS49.sim import step5_pad
step5_pad.big_data = ls49_big_data
from LS49.sim.step5_pad import full_path
# %%%%%%

if __name__=="__main__":

  from matplotlib import pyplot as plt

  GS = george_sherrell(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
  GS.plot_them(plt,f1="b.",f2="b.")
  GS.plot_them(plt,f1="b-",f2="b-")
  GS = george_sherrell(full_path("data_sherrell/pf-rd-red_fftkk.out"))
  GS.plot_them(plt,f1="r.",f2="r.")
  GS.plot_them(plt,f1="r-",f2="r-")
  GS = george_sherrell(full_path("data_sherrell/Fe.dat"))
  GS.plot_them(plt,f1="m-",f2="m-")

  from scipy.interpolate import interp1d
  fFe0fp = interp1d(GS.energy, GS.fp, kind="cubic")
  fFe0fdp= interp1d(GS.energy, GS.fdp, kind="linear")
  xnew = range(7070,7180)
  plt.plot(xnew, fFe0fp(xnew), "g+")
  plt.plot(xnew, fFe0fdp(xnew), "g+")
  for item in xnew:
    print ("%11.2f%15.7f%15.7f"%(float(item),fFe0fp(item),fFe0fdp(item)))

  plt.axes().set_xlim((7088,7152))
  plt.axes().set_ylim((-8.3,4.2))


  print(list(GS.energy))
  print(list(GS.fdp))

  plt.show()
