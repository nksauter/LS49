#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division, unicode_literals

from scitbx.array_family import flex

from LS49.work2_for_aca_lsq.abc_background import fit_roi_multichannel # implicit import
# multichannel needed for unpickling

# %%% boilerplate specialize to packaged big data %%%
import os
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
# %%%%%%
from LS49.sim.step5_pad import data
local_data = data()
Fe_oxidized_model = local_data.get("Fe_oxidized_model")
Fe_reduced_model = local_data.get("Fe_reduced_model")
Fe_metallic_model = local_data.get("Fe_metallic_model")

from LS49.ML_push.new_global_fdp_refinery import MPI_Run as upstream_base_script
from LS49.sim.fdp_plot import george_sherrell
from LS49.ML_push.new_global_fdp_refinery import rank_0_fit_all_f, george_sherrell_star

def plot_em(self,key,values):
    self.x = values # XXX
    if not self.plot_plt_imported:
      from matplotlib import pyplot as plt
      self.plt = plt
      if self.params.LLG_evaluator.title is None:
        self.plt.ion() # interactive - on
      self.plot_plt_imported = True
    if self.params.LLG_evaluator.title is None:
      self.plt.cla() #clear last access
    fine = self.params.LLG_evaluator.plot_interpolation # plot the non-modeled f values
    fig = self.plt.figure()

    # ground truth
    from LS49.sim.step5_pad import full_path
    GS = george_sherrell(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
    GS.plot_them(self.plt,f1="b-",f2="b-")
    GS = george_sherrell(full_path("data_sherrell/pf-rd-red_fftkk.out"))
    GS.plot_them(self.plt,f1="r-",f2="r-")
    GS = george_sherrell(full_path("data_sherrell/Fe_fake.dat")) # with interpolated points
    GS.plot_them(self.plt,f1="m-",f2="m-")

    # starting values
    GS = george_sherrell_star(fp = self.starting_params_FE1[0:100],fdp = self.starting_params_FE1[100:200])
    GS.plot_them(fine,self.plt,f1="bx",f2="bx")
    GS = george_sherrell_star(fp = self.starting_params_FE2[0:100],fdp = self.starting_params_FE2[100:200])
    GS.plot_them(fine,self.plt,f1="rx",f2="rx")

    # current values
    GS = george_sherrell_star(fp = self.x[0:100],fdp = self.x[100:200])
    GS.plot_them(fine,self.plt,f1="b.",f2="b.")
    GS = george_sherrell_star(fp = self.x[200:300],fdp = self.x[300:400])
    GS.plot_them(fine,self.plt,f1="r.",f2="r.")

    self.plt.axes().set_xlim((7102,7137)) # XXX 7088,7152
    self.plt.axes().set_ylim((-8.6,4.5))
    self.plt.title("Macrocycle %d Iteration %d"%(self.macrocycle,self.iteration)) # XXX
    if self.params.LLG_evaluator.title is not None:
      macrocycle_tell = "" if self.macrocycle is None else "macrocycle_%02d_"%self.macrocycle
      fig.savefig("replot_%s_%siteration_%02d.png"%(self.params.LLG_evaluator.title,
                  macrocycle_tell,self.iteration))
    else:
      self.plt.draw()
      self.plt.pause(0.2)
    fig.clf() # clear figure XXX
    #self.plt.show()

def plot_em_broken(self,key,values):
    self.x = values # XXX
    if not self.plot_plt_imported:
      from matplotlib import pyplot as plt
      self.plt = plt
      if self.params.LLG_evaluator.title is None:
        self.plt.ion() # interactive - on
      self.plot_plt_imported = True
    if self.params.LLG_evaluator.title is None:
      self.plt.cla() #clear last access
    fine = self.params.LLG_evaluator.plot_interpolation # plot the non-modeled f values
    fig, (ax1, ax2) = self.plt.subplots(2,1,sharex=True,squeeze=True)

    # ground truth
    from LS49.sim.step5_pad import full_path
    GS = george_sherrell(full_path("data_sherrell/pf-rd-ox_fftkk.out"))
    GS.plot_them(ax1,f1="b-",f2="b-")
    GS.plot_them(ax2,f1="b-",f2="b-")
    GS = george_sherrell(full_path("data_sherrell/pf-rd-red_fftkk.out"))
    GS.plot_them(ax1,f1="r-",f2="r-")
    GS.plot_them(ax2,f1="r-",f2="r-")
    GS = george_sherrell(full_path("data_sherrell/Fe_fake.dat")) # with interpolated points
    GS.plot_them(ax1,f1="m-",f2="m-")
    GS.plot_them(ax2,f1="m-",f2="m-")

    # starting values
    GS = george_sherrell_star(fp = self.starting_params_FE1[0:100],fdp = self.starting_params_FE1[100:200])
    GS.plot_them(fine,ax1,f1="bx",f2="bx")
    GS.plot_them(fine,ax2,f1="bx",f2="bx")
    GS = george_sherrell_star(fp = self.starting_params_FE2[0:100],fdp = self.starting_params_FE2[100:200])
    GS.plot_them(fine,ax1,f1="rx",f2="rx")
    GS.plot_them(fine,ax2,f1="rx",f2="rx")

    # current values
    GS = george_sherrell_star(fp = self.x[0:100],fdp = self.x[100:200])
    GS.plot_them(fine,ax1,f1="b.",f2="b.")
    GS.plot_them(fine,ax2,f1="b.",f2="b.")
    GS = george_sherrell_star(fp = self.x[200:300],fdp = self.x[300:400])
    GS.plot_them(fine,ax1,f1="r.",f2="r.")
    GS.plot_them(fine,ax2,f1="r.",f2="r.")

    #self.plt.axes().set_xlim((7102,7137)) # XXX 7088,7152
    ax1.set_xlim(7102,7138)
    ax2.set_xlabel("Energy (eV)")
    ax1.set_ylabel("∆ f ′′")
    ax2.set_ylabel("∆ f ′")
    ax2.set_ylim(-8.6,-5.1)
    ax1.set_ylim(0.1, 4.5)
    ax1.set_title("Macrocycle %d Iteration %d"%(self.macrocycle,self.iteration)) # XXX
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False,linewidth=1)
    ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    if self.params.LLG_evaluator.title is not None:
      macrocycle_tell = "" if self.macrocycle is None else "macrocycle_%02d_"%self.macrocycle
      fig.savefig("replot_%s_%siteration_%02d.png"%(self.params.LLG_evaluator.title,
                  macrocycle_tell,self.iteration))
    else:
      self.plt.draw()
      self.plt.pause(0.2)
    #fig.clf() # clear figure XXX
    #self.plt.show()
rank_0_fit_all_f.plot_em = plot_em_broken


class MPI_Run(upstream_base_script):
  def data_from_log(self,logpath):
    with open(logpath,"r") as F:
      results = {}
      for line in F:
        if line.find("Macrocycle")==0:
          tokens = line.strip().split()
          key = int(tokens[1]),int(tokens[3]) # macrocycle, iteration
          values_str = [t.replace(",","").replace("[","").replace("]","") for t in tokens[4:]]
          values_float = flex.double([float(t) for t in values_str])
          results[key] = values_float
    return results


  def run(self):
    self.parse_input()
    W = rank_0_fit_all_f( self.params,
                          FE1_model=local_data.get(self.params.starting_model.preset.FE1),
                          FE2_model=local_data.get(self.params.starting_model.preset.FE2))
    W.reinitialize(logical_rank=0, comm_size=1, per_rank_items=None, per_rank_keys=None, per_rank_G=None,
                   HKL_lookup=None, static_fcalcs=None, model_intensities=None)

    logpath = "%s.out"%self.params.LLG_evaluator.title
    data = self.data_from_log(logpath)
    W.plot_plt_imported = False
    for ky,v in data.iteritems():
      W.iteration = ky[1]
      W.macrocycle = ky[0]
      W.plot_em(ky,v)

    print(logpath,"ended at iteration")
    self.mpi_helper.comm.barrier()

if __name__=="__main__":
  Usage = """
libtbx.python ../modules/LS49/ML_push/replot.py LLG_evaluator.enable_plot=True starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model  LLG_evaluator.title=metal_metal_16k LLG_evaluator.plot_interpolation=True

libtbx.python ../modules/LS49/ML_push/replot.py LLG_evaluator.enable_plot=True LLG_evaluator.title=metal_metal_16k
convert -delay 12 -loop 0 *.png Fe_metal_iteration.gif
  """

  script = MPI_Run()
  result = script.run()
