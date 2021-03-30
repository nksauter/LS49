from __future__ import division
from dxtbx.model.experiment_list import ExperimentListFactory
import matplotlib.pyplot as plt
import time
import numpy as np

def binning_detail(x_input,y_input):
  minx = int(x_input[0]) + 3.5
  maxx = int(x_input[-1]) - 3.5

  bases = {}
  counts = {}
  heights = {}
  for idx in range(len(x_input)):
    xenergy = x_input[idx]
    yintensity = y_input[idx]
    stdbin = int(xenergy) + 0.5
    if stdbin < minx: continue
    if stdbin > maxx: continue
    bases[stdbin] = bases.get(stdbin,0.) + (x_input[idx+1] - x_input[idx-1])/2.0
    counts[stdbin] = counts.get(stdbin, 0.) + 1.0
    heights[stdbin] = heights.get(stdbin,0.) + yintensity
  x_result = np.array(list(bases.keys()))
  mean_count = np.array(list(counts.values())).mean()
  y_result = np.array([bases[k]*heights[k]/mean_count for k in bases])
  # finally kludge the intensities so they never drop below zero
  bottom = min(y_result)
  if bottom < 0.0: y_result -= bottom
  return x_result, y_result

def method3(energies_raw, weights_raw):
    lower = weights_raw[0:50].mean()
    upper = weights_raw[-50:].mean()
    baseline = (np.array(range(len(weights_raw)))/len(weights_raw))*(upper-lower)+lower
    real = np.array(list(weights_raw-baseline))
    #from IPython import embed; embed()
    #plt.plot(xrange(len(real)), real, 'r-')
    #plt.show()
    fr = np.fft.rfft(real)
    #print(type(fr), len(real)//2, len(fr.real))
    #plt.plot(xrange(len(fr.real)), fr.real, 'b-')
    #plt.show()
    #low_pass_fr:
    for x in range(len(fr)//4, len(fr)):
      fr[x]=0.+0.j
    filtered_real = np.fft.irfft(fr)

    fit_x,fit_y = binning_detail(energies_raw, filtered_real)
    return fit_x, fit_y, filtered_real

if __name__=="__main__":
  exptfile = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split_c/split_%04d.expt"
  init=True
  for idx in range(600,700):
    experiment_file = exptfile%idx

    El = ExperimentListFactory.from_json_file(experiment_file,
                                              check_format=True)
    exper = El[0]
    beam = exper.beam
    spec = exper.imageset.get_spectrum(0)
    energies_raw, weights_raw = spec.get_energies_eV().as_numpy_array(), \
                                spec.get_weights().as_numpy_array()

    fit_x, fit_y, filtered_real = method3(energies_raw, weights_raw)
    if init:
      plt.ion()
      figure, ax = plt.subplots(figsize=(8,6))
      line1, = ax.plot(energies_raw, weights_raw, "b.")
      line2, = ax.plot(energies_raw, filtered_real, "r-")
      line3, = ax.plot(fit_x, fit_y, "g.")
      plt.title("SwissFEL X-ray pulses",fontsize=25)
      plt.xlabel("X-ray energy (eV)",fontsize=18)
      plt.ylabel("intensity",fontsize=18)
      plt.ylim(0,30000)
      init=False
    if not init:
      updated_y = weights_raw
      line1.set_xdata(energies_raw)
      line1.set_ydata(weights_raw)
      line2.set_xdata(energies_raw)
      line2.set_ydata(filtered_real)
      line3.set_xdata(fit_x)
      line3.set_ydata(fit_y)
      figure.canvas.draw()
      figure.canvas.flush_events()
      #input()
      #time.sleep(0.1)
