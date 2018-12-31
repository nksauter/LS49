from __future__ import print_function, division
from six.moves import range
from six.moves import cPickle as pickle
from scitbx.matrix import col,sqr
from dials.algorithms.shoebox import MaskCode
from scitbx.array_family import flex
import numpy as np
import math
from matplotlib import pyplot as plt

json_glob = "LS49_integ_step5/idx-step5_MPIbatch_0%05d.img_integrated_experiments.json"
pickle_glob = "LS49_integ_step5/idx-step5_MPIbatch_0%05d.img_integrated.pickle"

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

#specialize this file to look at one particular index
distance_mm = 141.7
pixel_sz_mm = 0.11
mos_rotation_deg = 0.05

plot = False

def plot_energy_scale_noplot(SS,d_Ang,abs_PA,origin,position0,B,intensity_lookup,intensity_lookup_1,key):
  unit_pos0 = position0.normalize()
  spectrumx = []
  spectrumy = []
  spectrumy_1 = flex.double()
  for eV in range(7090,7151):
    spectrumx.append(eV)
    specy = 0.
    specy_1 = 0.
    lambda_Ang = 12398.425 / eV
    two_theta = 2. * math.asin( lambda_Ang / (2.*d_Ang))
    radius_mm = distance_mm * math.tan(two_theta)
    radius_px = radius_mm / pixel_sz_mm

    for rot in range(-8,9):
      PA = abs_PA + 0.25*rot*mos_rotation_deg*math.pi/180.
      clock = unit_pos0.rotate_2d(-PA, deg=False)
      position1 = origin + radius_px*clock
      int_coords = (int(position1[0]),int(position1[1]))
      specy += intensity_lookup.get(int_coords,0)
      specy_1 += intensity_lookup_1.get(int_coords,0)
    spectrumy.append(specy)
    spectrumy_1.append(specy_1)

  iterator = SS.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)
  wavlen, flux, wavelength_A = next(iterator) # list of lambdas, list of fluxes, average wavelength
  # the "flux" is the spectrum from the FEE spectrometer (simulated)

  combined_model = flex.double()
  incident_xaxis = 12398.425/wavlen
  int_ix = [int (ix) for ix in incident_xaxis]
  for ic in range(len(spectrumx)):
    ic_idx = int_ix.index(spectrumx[ic])
    combined_model.append(flux[ic_idx] * spectrumy_1[ic])
  cscale = max(spectrumy)/max(combined_model)
  CC=flex.linear_correlation(combined_model, flex.double(spectrumy)).coefficient()
  print ("The correlation coefficient is",CC)
  return spectrumx,spectrumy,combined_model,CC
  # spectrumy is the observed spectrum, projected
  # spectrumy_1 is the partiality model x Icalc, projected
  # combined model is the partiality x Icalc x spectrum, projected

def parse_postrefine():
  #lines = open("/net/dials/raid1/sauter/LS49_merge/merge5_redo2.log")
  lines = open("./LS49_merge/merge5_redo2.log")
  result = {}
  for line in lines:
    if "ASTAR" not in line: continue
    stripped = line.replace("(","").replace(")","").replace(",","")
    tokens = stripped.split()
    try:
      serial_no = int(tokens[1][-17:-11])
      rotmat = sqr(tuple([float(a) for a in tokens[2:]]))
      assert len(rotmat)==9
      result[serial_no] = rotmat
    except (IndexError,ValueError,AssertionError):
      continue
  return result
  # a lot of files have dropped out because of thread-to-thread contention in the log
  # file (thousands).   Should re-run this sometime with nproc=1

def get_items(myrank=None,mykey=None):
  postreffed = parse_postrefine()
  print ("# postrefined images",len(postreffed))
  maxy = None
  ycount = 0
  for key in postreffed:
    if mykey is not None:
      if key != mykey:continue
    #each rank should only allow keys in the range from myrank*N_stride to (myrank+1)*N_stride
    if key<myrank*N_stride: continue
    if key >= (myrank+1) * N_stride: continue
    if key >= N_total : continue

    from dxtbx.model.experiment_list import ExperimentListFactory
    E = ExperimentListFactory.from_json_file(json_glob%key,check_format=False)[0]
    C = E.crystal
    C.show()
    from six.moves import cPickle as pickle
    T = pickle.load(open(pickle_glob%key,"rb"))
    resolutions = T["d"]
    millers = T["miller_index"]
    nitem = len(resolutions)
    ycount+=1
    if maxy is not None and ycount>maxy:
      print("maxy break %d"%maxy)
      break
    print ("THE ACTUAL JSON / PICKLE USED:",json_glob%key,pickle_glob%key)
    yield T,key

class fit_background_abc_lsq(object):
  def __init__(self,shoebox,verbose=False):
    self.sb = shoebox
    data = self.sb.data - self.sb.background
    F = sb.data.focus()
    sum_x_sq = 0.; sum_x = 0.; sum_y_sq = 0.; sum_y = 0.;sum_xy = 0.; sum_xd = 0.;sum_yd = 0.; sum_d = 0.;N=0
    for x in range(F[1]):
      for y in range(F[2]):
        if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
          pass # not part of background
        else:
          # could calculate analytically, all nine elements? No, not with mask
          sum_x_sq += x*x; sum_y_sq += y*y; sum_x += x; sum_y += y; sum_xy += x*y ; N += 1
          datapt = self.sb.data[0,x,y]
          sum_xd += x * datapt; sum_yd += y*datapt; sum_d += datapt
    # Now solve equation.
    LSQ_matrix = sqr((sum_x_sq,sum_xy,sum_x,sum_xy,sum_y_sq,sum_y,sum_x,sum_y,N))
    RHS = col((sum_xd, sum_yd, sum_d))
    ABC = LSQ_matrix.inverse() * RHS
    self.abc = ABC

    if verbose:
      # Background subtracted data
      for x in range(F[1]):
        for y in range(F[2]):
          print ("%4.0f"%(data[0,x,y]),end=' ')
        print()
      print()

      # Same, with DIALS foreground mask shown
      for x in range(F[1]):
        for y in range(F[2]):
          if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
            print ("   M",end=' ')
          elif sb.mask[c]&MaskCode.Valid != MaskCode.Valid:
            print ("   V",end=' ')
          else:
            print ("%4.0f"%(data[0,x,y]),end=' ')
        print()
      print()

      # Raw data, with DIALS foreground mask shown
      for x in range(F[1]):
        for y in range(F[2]):
          if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
            print ("   M",end=' ')
          else:
            print ("%4.0f"%(self.sb.data[0,x,y]),end=' ')
        print()
      print()

      # LSQ background model
      for x in range(F[1]):
        for y in range(F[2]):
            print ("%4.0f"%(ABC[0]*x+ABC[1]*y+ABC[2]),end=' ')
        print()
      print()

      # DIALS background model
      for x in range(F[1]):
        for y in range(F[2]):
          print ("%4.0f"%(self.sb.background[0,x,y]),end=' ')
        print()
      print()
    # analyse the variance
    print (self.abc.elems)
    diff = flex.double()
    approx_poisson_sigma = flex.double()
    for x in range(F[1]):
      for y in range(F[2]):
        if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
          pass # not part of background
        else:
          diff.append(data[0,x,y])
          approx_poisson_sigma.append(math.sqrt(ABC[0]*x+ABC[1]*y+ABC[2]))
    print ("Min=%4.0f, Max=%4.0f"%(flex.min(diff), flex.max(diff)) )
    MV = flex.mean_and_variance(diff)
    print ("mean=%4.0f, stddev=%4.0f"%(MV.mean(),MV.unweighted_sample_standard_deviation()))
    print ("mean Poisson stddev=%4.0f"%(flex.mean(approx_poisson_sigma)))

class fit_background_abc_ml(object):
  def __init__(self,shoebox,verbose=False,basic_verbose=False):
    import scitbx
    self.sb = shoebox
    self.x = flex.double([0.,0.,100.,]) # setting background model paramteres to initial guess zero
    self.n = 3 # three parameters
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.a = self.x
    diffdata = self.sb.data - self.sb.background
    F = self.sb.data.focus()

    if basic_verbose:
      # Background subtracted data
      for x in range(F[1]):
        for y in range(F[2]):
          print ("%4.0f"%(diffdata[0,x,y]),end=' ')
        print()
      print()

      # Same, with DIALS foreground mask shown
      for x in range(F[1]):
        for y in range(F[2]):
          if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
            print ("   M",end=' ')
          else:
            print ("%4.0f"%(diffdata[0,x,y]),end=' ')
        print()
      print()

      # Raw data, with DIALS foreground mask shown
      for x in range(F[1]):
        for y in range(F[2]):
          if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
            print ("   M",end=' ')
          else:
            print ("%4.0f"%(self.sb.data[0,x,y]),end=' ')
        print()
      print()

      # DIALS background model
      for x in range(F[1]):
        for y in range(F[2]):
          print ("%4.0f"%(self.sb.background[0,x,y]),end=' ')
        print()
      print()

    if verbose:
      # ML background model
      for x in range(F[1]):
        for y in range(F[2]):
            print ("%4.0f"%(self.a[0]*x+self.a[1]*y+self.a[2]),end=' ')
        print()
      print()

    # analyse the variance
    print (list(self.a))
    diff = flex.double()
    approx_poisson_sigma = flex.double()
    for x in range(F[1]):
      for y in range(F[2]):
        if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
          pass # not part of background
        else:
          diff.append(diffdata[0,x,y])
          approx_poisson_sigma.append(math.sqrt(self.a[0]*x+self.a[1]*y+self.a[2]))
    print ("ML Min=%4.0f, Max=%4.0f"%(flex.min(diff), flex.max(diff)) )
    MV = flex.mean_and_variance(diff)
    print ("ML mean=%4.0f, stddev=%4.0f"%(MV.mean(),MV.unweighted_sample_standard_deviation()))
    print ("ML mean Poisson stddev=%4.0f"%(flex.mean(approx_poisson_sigma)))

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%10.4f"%a for a in self.x]),"]")

  def compute_functional_and_gradients(self):
    self.a = self.x
    f = 0.; ga = 0.; gb = 0.; gc = 0.
    F = self.sb.data.focus()
    for x in range(F[1]):
      for y in range(F[2]):
        if self.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
          pass # not part of background
        else:
          model_lambda = self.a[0]*x+self.a[1]*y+self.a[2]
          if model_lambda<=0:
            f+= model_lambda # complete kludge, guard against math domain error
          else:
            datapt = self.sb.data[0,x,y]
            f += model_lambda - datapt * math.log(model_lambda)
          ga += x * (1. - datapt/model_lambda) # from handwritten notes
          gb += y * (1. - datapt/model_lambda)
          gc += (1. - datapt/model_lambda)
    self.print_step("LBFGS stp",f)
    return f, flex.double([ga,gb,gc])

class fit_roi(object):
  def __init__(self,bkgrd,roi,verbose=True,basic_verbose=True):
    import scitbx
    self.sb = bkgrd.sb # shoebox
    self.roi = roi
    self.x = flex.double([bkgrd.a[0],bkgrd.a[1],bkgrd.a[2],1.]
             ) # setting background model paramteres to previously determined values
    self.n = 4 # four parameters
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.a = self.x
    diffdata = self.sb.data - self.sb.background
    F = self.sb.data.focus()

    if basic_verbose:
      # Background subtracted data
      for x in range(F[1]):
        for y in range(F[2]):
          print ("%4.0f"%(diffdata[0,x,y]),end=' ')
        print()
      print()

      # Same, with DIALS foreground mask shown
      for x in range(F[1]):
        for y in range(F[2]):
          if bkgrd.sb.mask[0,x,y]&MaskCode.Foreground == MaskCode.Foreground:
            print ("   M",end=' ')
          else:
            print ("%4.0f"%(diffdata[0,x,y]),end=' ')
        print()
      print()

      # Raw data, with DIALS foreground mask shown
      for x in range(F[1]):
        for y in range(F[2]):
            print ("%4.0f"%(bkgrd.sb.data[0,x,y]),end=' ')
        print()
      print()

    if verbose:
      # ML background + ML roi model
      for x in range(F[1]):
        for y in range(F[2]):
          background = self.a[0]*x+self.a[1]*y+self.a[2]
          print ("%4.0f"%(background + self.a[3] * self.roi[x,y]),end=' ')
        print()
      print()

      # bkgrd background model (before modeling diffraction)
      for x in range(F[1]):
        for y in range(F[2]):
          background = bkgrd.a[0]*x+bkgrd.a[1]*y+bkgrd.a[2]
          print ("%4.0f"%(background),end=' ')
        print()
      print()

      print ("# Modelled ROI")
      for x in range(F[1]):
        for y in range(F[2]):
            print ("%4.0f"%(self.a[3] * roi[x,y]),end=' ')
        print()
      print()

      self.show_residual(plot=False)

  def show_residual(self,plot=False):
    F = self.sb.data.focus()
    if plot:
      # Residual
      for x in range(F[1]):
        for y in range(F[2]):
          background = self.a[0]*x+self.a[1]*y+self.a[2]
          model = background + self.a[3] * self.roi[x,y]
          print ("%4.0f"%(self.sb.data[0,x,y]-model),end=' ')
        print()
      print()

    # analyse the variance
    print (list(self.a))
    diff = flex.double()
    approx_poisson_sigma = flex.double()
    for x in range(F[1]):
      for y in range(F[2]):
          background = self.a[0]*x+self.a[1]*y+self.a[2]
          model = background + self.a[3] * self.roi[x,y]
          diffdata = self.sb.data[0,x,y]-model
          diff.append(diffdata)
          approx_poisson_sigma.append(math.sqrt(model))
    MV = flex.mean_and_variance(diff)
    fmin,fmax = flex.min(diff), flex.max(diff)
    print ("residual Min=%4.0f, Mean=%4.0f, Max=%4.0f"%(fmin, MV.mean(), fmax) )
    print ("residual stddev=%4.0f"%(MV.unweighted_sample_standard_deviation()))
    print ("ML mean Poisson stddev=%4.0f"%(flex.mean(approx_poisson_sigma)))
    if plot:
      from matplotlib import pyplot as plt
      n, bins, patches = plt.hist(diff, 40 ,normed=1, facecolor='g', alpha=0.75)
      plt.xlabel('Model vs. Observation Residual')
      plt.ylabel('Probability')
      plt.title('Histogram of Model Residual')
      plt.show()


  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%10.4f"%a for a in self.x]),"]")

  def compute_functional_and_gradients(self):
    self.a = self.x
    f = 0.; ga = 0.; gb = 0.; gc = 0. ; gG = 0.
    F = self.sb.data.focus()
    for x in range(F[1]):
      for y in range(F[2]):
          model_lambda = self.a[0]*x+self.a[1]*y+self.a[2]+self.a[3]*self.roi[x,y]
          if model_lambda<=0:
            f+= model_lambda # complete kludge, guard against math domain error
          else:
            datapt = self.sb.data[0,x,y]
            f += model_lambda - datapt * math.log(model_lambda)
          ga += x * (1. - datapt/model_lambda) # from handwritten notes
          gb += y * (1. - datapt/model_lambda)
          gc += (1. - datapt/model_lambda)
          gG += self.roi[x,y] * (1. - datapt/model_lambda)
    self.print_step("LBFGS stp",f)
    return f, flex.double([ga,gb,gc,gG])

class fit_roi_multichannel(object):
  def __init__(self,bkgrd,channels,rescale,verbose=True,basic_verbose=True):
    import scitbx
    self.sb_data = bkgrd.sb.data # shoebox
    self.bkgrd_a = bkgrd.a
    self.channels = channels
    self.roi = rescale[0] * self.channels[0]
    assert (len(self.channels) == 50) # there are only 50 channels, spaced 2 eV
    for ichannel in range(1,len(self.channels)):
      if flex.sum(self.channels[2*ichannel])>50:
        print (ichannel, 2*ichannel, rescale[2*ichannel])
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

from LS49.work2_for_aca_lsq.util_partiality import get_partiality_response

if __name__=="__main__":

  Usage = """mpirun -n 50 libtbx.python gen_data_mpi.py
             rather: for x in `seq 0 64`; do libtbx.python jun24_gen_data_mpi.py $x & done
             for x in `seq 0 3`; do time libtbx.python may27_gen_data_mpi.py $x > /dev/null & done"""
  from libtbx.mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()
  N_total = 100000 # number of items to simulate
  usingMPI = True
  if rank==0 and size==1: # special case of testing it
    import sys
    try:
      rank = int (sys.argv[1])
      size=1024
      usingMPI = False
    except Exception: pass
    print("MPI rank %d of %d"%(rank,size))
  N_stride = int(math.ceil(N_total/size)) # total number of tasks per rank
  print ("hello from rank %d of %d with stride %d"%(rank,size,N_stride))

  if (not usingMPI) or rank == 0:
    print ("set up in rank 0")
    from LS49.sim.step5_pad import data
    pdb_lines = data().get("pdb_lines")
    from LS49.sim.util_fmodel import gen_fmodel

    GF = gen_fmodel(resolution=10.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.7)
    A = GF.get_amplitudes()

    from LS49.spectra.generate_spectra import spectra_simulation
    SS = spectra_simulation()

    ori_N_total = N_total # number of items to simulate
    mt = flex.mersenne_twister(seed=0)
    random_orientations = []
    for iteration in range(ori_N_total):
      random_orientations.append( mt.random_double_r3_rotation_matrix() )

    transmitted_info = dict(spectra = SS,
                            amplitudes = A,
                            orientations = random_orientations,
    )
    with (open("confirm_P1_range_reduced_intensities_dict.pickle","rb")) as F: # Einsle reduced
    #with (open("confirm_P1_range_oxidized_intensities_dict.pickle","rb")) as F: # Einsle oxidized
    #with (open("confirm_P1_range_metallic_intensities_dict.pickle","rb")) as F: # Einsle metallic
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

  origin = col((1500,1500))
  position0 = col((1500,3000))-origin
  nitem = 0
  nall_spots = 0
  nres_range = 0
  npos_angle = 0
  nVF = 0
  millerd = {}
  intensity_dict = transmitted_info["intensity_dict"]
  #for item,key in get_items(key=3271):
  for item,key in get_items(rank):
    result = dict(image=key,asu_idx_C2_setting=[],orig_idx_C2_setting=[],pick=pickle_glob%key,spectrumx=[],obs=[],model=[],cc=[],
                  simtbx_millers=[],simtbx_intensity=[])
    roi_results = []
    #good indices: 3271, 2301: continue
    d = item["d"]
    nitem += 1

    print ("ENTER key %d in rank %d"%(key,rank))

    nall_spots += len(item)
    iselect = ((d < 2.5) & (d > 2.1))
    nres_range += len(d.select(iselect))

    # geometric selection:  between position angles 150 and 210 degrees.
    hkl = item["miller_index"].select(iselect)
    cust_copy = transmitted_info["amplitudes"].customized_copy(indices=hkl,data=flex.double(len(hkl)),sigmas=flex.double(len(hkl)))
    asu = cust_copy.map_to_asu().indices()

    # the indices are already in the 2.1 to 2.5 range
    # hkl is original index; asu is asu index

    xyz = item["xyzobs.px.value"].select(iselect)
    calcpx = item["xyzcal.px"].select(iselect)
    shoe = item["shoebox"].select(iselect)
    intensity_lookup ={}
    intensity_lookup_1 ={}

    for x in range(len(hkl)):
      slow = xyz[x][1]
      fast = xyz[x][0]
      positionX = col((slow,fast))-origin
      position_angle = positionX.angle(position0,deg=True)
      if position_angle > 150.:
        print ("PA key %d in rank %d hkl %s"%(key,rank,hkl[x]))

        npos_angle += 1
        millerd[asu[x]]=millerd.get(asu[x],0)+1
        sb = shoe[x]
        nsb = sb.mask.size()
        for c in range(nsb):
          if sb.mask[c]&MaskCode.Valid == MaskCode.Valid and sb.mask[c]&MaskCode.Foreground == MaskCode.Foreground:
            nVF += 1
          intensity_lookup[(int(sb.coords()[c][1]),int(sb.coords()[c][0]))] = sb.data[c]-sb.background[c]
        spotprediction = calcpx[x] # DIALS coords (fast,slow)
        spotvec = col((spotprediction[0],spotprediction[1]))-origin # panel coords (center=origin)
        abs_PA = math.atan2(spotvec[1],spotvec[0]) # clockwise plotted on image_viewer (because vertical axis upside down)
        B = sb.bbox
        ROI = ((B[0],B[1]),(B[2],B[3]))
        try:
          fb_lsq = fit_background_abc_lsq(sb,verbose=False)
          fb_ml  = fit_background_abc_ml(sb,verbose=False)
        except Exception as e:
          print ("FAILing fit_background_abc on",e)
          continue # avoid LBFGS Runtime and other exceptions
        values = sb.data-sb.background # ADU above background
        # intensity_lookup consists of the "observed" data from shoeboxes

        v0 = values.set_selected(values<=0, 0.)
        v1 = v0.set_selected(v0>255,255)
        v2 = (256.-v1)/256.
        np_v2 = np.ndarray(shape=(B[3]-B[2],B[1]-B[0],), dtype=np.float32, buffer=v2.as_numpy_array())

        # insert code here to estimate the partiality response
        PRD = get_partiality_response(key,hkl[x],spectra_simulation=transmitted_info["spectra"],ROI=ROI,
                                      rand_ori = sqr(transmitted_info["orientations"][key]),
                                      tophat_spectrum=False,
                                     )
        pr_value=PRD["roi_pixels"]
        miller=PRD["miller"]
        intensity=PRD["intensity"]
        channels = PRD["channels"]

        FR = fit_roi(fb_ml, pr_value)

        # commented out code proves that pr_value is the sum of the channels
        if False:
          F = sb.data.focus()
          for x in range(F[1]):
            for y in range(F[2]):
              print ("%4.0f"%(FR.a[3] * pr_value[x,y]),end=' ')
            print()
          print()
          monoband_residual = pr_value.deep_copy()
          for ckey in channels:
            monoband_residual -= channels[ckey]
            print("monoband residual %d"%ckey)
            for x in range(F[1]):
              for y in range(F[2]):
                print ("%4.0f"%(FR.a[3] * monoband_residual[x,y]),end=' ')
              print()
            print()
        if intensity==0.:
          print ("FAILing on divide by zero")
          continue
        rescale_factor = intensity_dict[miller]/intensity
        if False:
          from matplotlib import pyplot as plt
          plt.plot(range (len(rescale_factor)),rescale_factor,'r-')
          plt.show()
        try:
          FRC = fit_roi_multichannel(fb_ml, channels, rescale_factor)
        except Exception as e:
          print ("FAILing fit_roi_multichannel on",e)
          continue
        # condition the FRC for pickling:
        if True: # pickle it
           del FRC.minimizer
           FRC.simtbx_intensity_7122 = intensity
           FRC.simtbx_P1_miller = miller
           FRC.orig_idx_C2_setting=hkl[x]
           FRC.asu_idx_C2_setting=asu[x]
           FRC.image_no = key
           roi_results.append(FRC)
        print ("FR.a",list(FR.a))
        print ("FRC.a",list(FRC.a))
        print ("# Modelled ROI with rescaling")
        F = sb.data.focus()
        for ix in range(F[1]):
          for iy in range(F[2]):
              print ("%4.0f"%(FRC.a[3] * FRC.roi[ix,iy] - FR.a[3]*FR.roi[ix,iy]),end=' ')
              #print ("%4.0f"%(FRC.roi[ix,iy] - FR.roi[ix,iy]),end=' ')
          print()
        print()


        print (""" *************************
 *********************
 *********************
 LLG Image %06d Miller %s NLL constant F = %9.1f     channels F = %9.1f
 *********************
 *********************
"""%(key, hkl[x], FR.compute_functional_and_gradients()[0], FRC.compute_functional_and_gradients()[0]))
        if False: # confirm the consistency of structure factor models
          print ("P1 miller index, presumably",miller,"intensity",intensity)
          print (list(intensity_dict[miller]))
          x_axis = [7070.5 + i for i in range(100)]
          from matplotlib import pyplot as plt
          plt.plot(x_axis,intensity_dict[miller],'r-')
          plt.plot([7122,],intensity,"b.")
          plt.show()

        for c in range(nsb):
          intensity_lookup_1[(int(sb.coords()[c][1]),int(sb.coords()[c][0]))] = pr_value[c]
        if len(intensity_lookup_1) != len(intensity_lookup):
          print ("FAILing on intensity lookup unequal", len(intensity_lookup_1),len(intensity_lookup))
          continue
        assert len(pr_value) == len(sb.data)
        # intensity_lookup_1 consists of partiality model data from posthoc simulator (partiality x Icalc)

        values_1 = pr_value # sb.data # ADU
        v0_1 = values_1.set_selected(values_1<=0, 0.)
        v1_1 = v0_1.set_selected(v0_1>255,255)
        v2_1 = (256.-v1_1)/256.
        np_v2_1 = np.ndarray(shape=(B[3]-B[2],B[1]-B[0],), dtype=np.float64, buffer=v2_1.as_numpy_array())

        ssplot=False # XXXYYY
        if ssplot:
          ax=ax1=ax2=None
          fig = plt.figure(figsize=(8,6))
          ax = plt.subplot2grid ((2,2),(0,0))
          ax1 = plt.subplot2grid((2,2),(0,1))
          ax2 = plt.subplot2grid((2,2),(1,0),colspan=2)
          ax.imshow(np_v2, cmap=plt.cm.gray, interpolation="nearest")
          ax1.imshow(np_v2_1, cmap=plt.cm.gray, interpolation="nearest")
          plt.title("key %d in rank %d"%(key,rank))
          plt.show()

        spectrumx,spectrumy,combined_model,CC =\
        plot_energy_scale_noplot(transmitted_info["spectra"],
                                 transmitted_info["amplitudes"].unit_cell().d(hkl[x]),
                                 abs_PA,origin,position0,B,intensity_lookup,intensity_lookup_1,key)
        result["simtbx_millers"].append(miller)
        result["simtbx_intensity"].append(intensity)
        result["asu_idx_C2_setting"].append(asu[x])
        result["orig_idx_C2_setting"].append(hkl[x])
        result["spectrumx"].append(spectrumx) # spectrumx: energy channels in units of eV
        result["obs"].append(spectrumy) # observations, presumably
        result["model"].append(combined_model) # model presumably consisting of partiality x spectrum x Icalc
        result["cc"].append(CC)
        if False and CC>0.7:
          plt.plot(spectrumx,1E10*flex.double(spectrumy),"b-")
          plt.plot(spectrumx,combined_model,"k-")
          plt.show()

    #recast the ROI miller indices in the centered (hi-sym) space group that was used for DIALS integration
    CBOP = cust_copy.change_of_basis_op_to_primitive_setting()
    from cctbx.array_family import flex as cflex
    M = cflex.miller_index(result["simtbx_millers"])
    P = cust_copy.customized_copy(indices=M)
    Q = P.change_basis(CBOP.inverse())
    result["simtbx_millers_DIALS_setting"]=list(Q.indices())
    result["cb_op_simtbx_to_dials"] = CBOP.inverse()

    # At this point we optionally restrict the data
    #print ("pickling key %d in rank %d"%(key,rank),result)
    #pickle.dump(result,open("dataX%04d.pickle"%rank,"ab"))
    print ("pickling key %d in rank %d"%(key,rank),roi_results)
    with open("abc_coverage/abcX%06d.pickle"%key,"ab") as F:
      pickle.dump(roi_results,F, pickle.HIGHEST_PROTOCOL)
