from __future__ import division, print_function
from scitbx.array_family import flex
# %%% boilerplate specialize to packaged big data %%%
import os
from LS49.adse13_196 import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
from LS49 import ls49_big_data
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
from LS49.sim.util_fmodel import gen_fmodel
from LS49.adse13_196.step5_pad import data
# %%%%%%
from libtbx.development.timers import Profiler

def create_cpu_channels(utilize):
  wavelength_A = 1.74 # general ballpark X-ray wavelength in Angstroms
  wavlen = flex.double([12398.425/(7070.5 + w) for w in range(utilize)])
  direct_algo_res_limit = 1.7

  local_data = data() # later put this through broadcast

  GF = gen_fmodel(resolution=direct_algo_res_limit,
      pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=wavelength_A)
  GF.set_k_sol(0.435)
  GF.make_P1_primitive()

  # Generating sf for my wavelengths
  sfall_channels = {}
  for x in range(len(wavlen)):
    GF.reset_wavelength(wavlen[x])
    GF.reset_specific_at_wavelength(
      label_has="FE1",
      tables=local_data.get("Fe_oxidized_model"),newvalue=wavlen[x])
    GF.reset_specific_at_wavelength(
      label_has="FE2",
      tables=local_data.get("Fe_reduced_model"),newvalue=wavlen[x])
    sfall_channels[x]=GF.get_amplitudes()
    print("CPU channel",x)
  return sfall_channels

def create_gpu_channels(cpu_channels,utilize):
  from libtbx.mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  devices_per_node = int(os.environ.get("DEVICES_PER_NODE", 1))
  this_device = rank%devices_per_node

  from simtbx.gpu import gpu_energy_channels
  gpu_channels_singleton = gpu_energy_channels (
    deviceId = this_device )

  assert gpu_channels_singleton.get_deviceID()==this_device
  print ("QQQ to gpu %d channels"%gpu_channels_singleton.get_nchannels(),"rank",rank)

  if gpu_channels_singleton.get_nchannels() == 0: # if uninitialized
    P = Profiler("Initialize the channels singleton rank %d, device %d"%(rank,this_device))
    for x in range(len(cpu_channels)):
      print("starting with ",x)
      gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, cpu_channels[x].indices(), cpu_channels[x].data())
      print("Finished sending to gpu %d channels"%gpu_channels_singleton.get_nchannels())
    del P
    assert len(cpu_channels)==utilize

def create_gpu_channels_one_rank(cpu_channels,utilize):
  this_device = 0

  from simtbx.gpu import gpu_energy_channels
  gpu_channels_singleton = gpu_energy_channels (
    deviceId = this_device )

  assert gpu_channels_singleton.get_deviceID()==this_device
  print ("QQQ to gpu %d channels"%gpu_channels_singleton.get_nchannels(),"one rank")

  if gpu_channels_singleton.get_nchannels() == 0: # if uninitialized
    P = Profiler("Initialize the channels singleton rank None, device %d"%(this_device))
    for x in range(len(cpu_channels)):
      print("starting with ",x)
      gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, cpu_channels[x].indices(), cpu_channels[x].data())
      print("Finished sending to gpu %d channels"%gpu_channels_singleton.get_nchannels())
    del P
    assert len(cpu_channels)==utilize

if __name__=="__main__":
  utilize=10
  CPU = create_cpu_channels(utilize)
  create_gpu_channels_one_rank(CPU,utilize)
  create_gpu_channels(CPU,utilize)
  print ("OK")
