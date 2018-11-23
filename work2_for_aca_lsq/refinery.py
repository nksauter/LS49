from __future__ import print_function, division
from six.moves import range
from six.moves import cPickle as pickle
from scitbx.array_family import flex
import math

json_glob = "/global/cscratch1/sd/nksauter/proj-e/LS49_integ_step5/idx-step5_MPIbatch_0%05d.img_integrated_experiments.json"
image_glob = "/net/dials/raid1/sauter/LS49/step5_MPIbatch_0%05d.img.gz"

#specialize this file to look at one particular index
distance_mm = 141.7
pixel_sz_mm = 0.11
mos_rotation_deg = 0.05

def get_items(myrank):
  for key in range(N_total):
    #each rank should only allow keys in the range from myrank*N_stride to (myrank+1)*N_stride
    if key<myrank*N_stride: continue
    if key >= (myrank+1) * N_stride: continue
    if key >= N_total : continue
    try:
      with open("abc_coverage/abcX%06d.pickle"%key,"rb") as F:
        T = pickle.load(F)
    except IOError:
      print("No file abc_coverage/abcX%06d.pickle"%key)
      continue
    yield T,key


class fit_roi_multichannel_repeat:
  def __init__(self,bkgrd_a,sb_data,channels,rescale,verbose=True,basic_verbose=True):
    import scitbx
    self.sb_data = sb_data # shoebox
    self.bkgrd_a = bkgrd_a
    self.channels = channels
    self.roi = rescale[0] * self.channels[0]
    assert (len(self.channels) == 50) # there are only 50 channels, spaced 2 eV
    for ichannel in range(1,len(self.channels)):
      #if flex.sum(self.channels[2*ichannel])>50:
      #  print (ichannel, 2*ichannel, rescale[2*ichannel])
      self.roi += rescale[2*ichannel] * self.channels[2*ichannel]
    self.x = flex.double([self.bkgrd_a[0],self.bkgrd_a[1],self.bkgrd_a[2],1.]
             ) # setting background model paramteres to previously determined values
    self.n = 4 # four parameters
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.a = self.x

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%10.4f"%a for a in self.x]),"]")

  def compute_functional_and_gradients(self):
    self.a = self.x
    f = 0.; ga = 0.; gb = 0.; gc = 0. ; gG = 0.
    F = self.sb_data.focus()
    for x in range(F[1]):
      for y in range(F[2]):
        model_lambda = self.a[0]*x+self.a[1]*y+self.a[2]+self.a[3]*self.roi[x,y]
        if model_lambda<=0:
          f+= model_lambda # complete kludge, guard against math domain error
        else:
          datapt = self.sb_data[0,x,y]
          f += model_lambda - datapt * math.log(model_lambda)
        ga += x * (1. - datapt/model_lambda) # from handwritten notes
        gb += y * (1. - datapt/model_lambda)
        gc += (1. - datapt/model_lambda)
        gG += self.roi[x,y] * (1. - datapt/model_lambda)
    #self.print_step("LBFGS stp",f)
    return f, flex.double([ga,gb,gc,gG])

class fit_one_image_multispot:
  def __init__(self,list_of_images):
    import scitbx
    #lay out the parameters.
    self.n_spots = len(list_of_images)
    self.n = 3*self.n_spots + 1
    self.x = flex.double()
    for ispot in range(self.n_spots):
      self.x.append(list_of_images[ispot].bkgrd_a[0])
      self.x.append(list_of_images[ispot].bkgrd_a[1])
      self.x.append(list_of_images[ispot].bkgrd_a[2])
    self.x.append(1.)
    self.roi = []
    for ispot in range(self.n_spots):
      intensity = list_of_images[ispot].simtbx_intensity_7122
      rescale_factor = intensity_dict[list_of_images[ispot].simtbx_P1_miller]/intensity
      channels = list_of_images[ispot].channels
      self.roi.append(rescale_factor[0] * channels[0])
      for ichannel in range(1,len(channels)):
        self.roi[ispot] += rescale_factor[2*ichannel] * channels[2*ichannel]
    self.sb_data = []
    for ispot in range(self.n_spots):
      self.sb_data.append(list_of_images[ispot].sb_data)
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.a = self.x

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%10.4f"%a for a in self.x]),"]")

  def compute_functional_and_gradients(self):
    self.a = self.x
    f = 0.;
    g = flex.double(self.n)
    for ispot in range(self.n_spots):
      F = self.sb_data[ispot].focus()
      for x in range(F[1]):
        for y in range(F[2]):
          model_lambda = self.a[3*ispot+0]*x+self.a[3*ispot+1]*y+self.a[3*ispot+2]+ \
                         self.a[-1]*self.roi[ispot][x,y]
          if model_lambda<=0:
            f+= model_lambda # complete kludge, guard against math domain error
          else:
            datapt = self.sb_data[ispot][0,x,y]
            f += model_lambda - datapt * math.log(model_lambda)
          g[3*ispot+0] += x * (1. - datapt/model_lambda) # from handwritten notes
          g[3*ispot+1] += y * (1. - datapt/model_lambda)
          g[3*ispot+2] += (1. - datapt/model_lambda)
          g[-1] += self.roi[ispot][x,y] * (1. - datapt/model_lambda)
    #self.print_step("LBFGS stp",f)
    return f, g


if __name__=="__main__":

  Usage = """mpirun -n 50 libtbx.python gen_data_mpi.py
             rather: for x in `seq 0 64`; do libtbx.python jun24_gen_data_mpi.py $x & done
             for x in `seq 0 3`; do time libtbx.python may27_gen_data_mpi.py $x > /dev/null & done"""
  usingMPI = True # XXXYYY
  if usingMPI:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
  else:
    import sys
    rank = int (sys.argv[1])
    #size=64
    size=1024
  N_total = 100000 # number of items to simulate
  N_stride = int(math.ceil(N_total/size)) # total number of tasks per rank
  print ("hello from rank %d of %d with stride %d"%(rank,size,N_stride))

  if (not usingMPI) or rank == 0:
    print ("set up in rank 0")
    transmitted_info = dict(
    )
    #with (open("confirm_P1_range_intensities_dict.pickle","rb")) as F: # Einsle reduced
    #with (open("confirm_P1_range_oxidized_intensities_dict.pickle","rb")) as F: # Einsle oxidized
    with (open("confirm_P1_range_metallic_intensities_dict.pickle","rb")) as F: # Einsle metallic
      intensity_dict = pickle.load(F)
      transmitted_info["intensity_dict"] = intensity_dict
    print ("finished setup in rank 0")
  else:
    transmitted_info = None
  if usingMPI:
    transmitted_info = comm.bcast(transmitted_info, root = 0)
    comm.barrier()
    import os
    host = os.environ["HOST"]
    print ("barrier from rank %d of %d"%(rank,size),host)

  intensity_dict = transmitted_info["intensity_dict"]

  mode="per-image"
  if mode=="per-hkl":
   for item,key in get_items(rank):
    print ("key %d in rank %d"%(key,rank))
    for x in range(len(item)):
        roi = item[x]
        intensity = roi.simtbx_intensity_7122
        rescale_factor = intensity_dict[roi.simtbx_P1_miller]/intensity
        FRC = fit_roi_multichannel_repeat(roi.bkgrd_a,roi.sb_data, roi.channels, rescale_factor)
        # condition the FRC for pickling:
        print ("roi.a",list(roi.a))
        print ("FRC.a",list(FRC.a))
        print (""" *******************************************************************
 LLG Image %06d Miller %s NLL constant F = %9.1f     channels F = %9.1f
******************************************
"""%(key, roi.orig_idx_C2_setting, roi.compute_functional_and_gradients()[0], FRC.compute_functional_and_gradients()[0]))

  elif mode=="per-image":
    for item,key in get_items(rank):
      print ("key %d in rank %d"%(key,rank))
      FOI = fit_one_image_multispot(list_of_images=item)
      print ("result",list(FOI.a))
      print (""" *******************************************************************
 LLG Image %06d on %d Bragg spots NLL    channels F = %9.1f
******************************************
"""%(key, len(item), FOI.compute_functional_and_gradients()[0]))
