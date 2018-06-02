from __future__ import print_function, division
from six.moves import range
from scitbx.matrix import col,sqr
from dials.algorithms.shoebox import MaskCode
from scitbx.array_family import flex
import numpy as np
import math
from matplotlib import pyplot as plt

json_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/idx-step5_MPIbatch_0%05d.img_integrated_experiments.json"
image_glob = "/net/dials/raid1/sauter/LS49/step5_MPIbatch_0%05d.img.gz"
pickle_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/idx-step5_MPIbatch_0%05d.img_integrated.pickle"

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
  lines = open("/net/dials/raid1/sauter/LS49_merge/merge5_redo2.log")
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

def get_items(myrank):
  postreffed = parse_postrefine()
  print ("# postrefined images",len(postreffed))
  maxy = 2001
  ycount = 0
  for key in postreffed:
    #each rank should only allow keys in the range from myrank*N_stride to (myrank+1)*N_stride
    if key<myrank*N_stride: continue
    if key >= (myrank+1) * N_stride: continue

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
    if ycount>maxy:
      print("maxy break %d"%maxy)
      break
    print ("THE ACTUAL JSON / PICKLE USED:",json_glob%key,pickle_glob%key)
    yield T,key

from LS49.sim.util_partiality import get_partiality_response

if __name__=="__main__":

  Usage = """mpirun -n 50 libtbx.python gen_data_mpi.py
             rather: for x in `seq 0 49`; do libtbx.python may27_gen_data_mpi.py $x & done
             for x in `seq 0 3`; do time libtbx.python may27_gen_data_mpi.py $x > /dev/null & done"""
  usingMPI = True
  if usingMPI:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
  else:
    import sys
    rank = int (sys.argv[1])
    size=50
  N_total = 100000 # number of items to simulate
  N_stride = 2000 # total number of tasks per rank
  print ("hello from rank %d of %d"%(rank,size))

  if (not usingMPI) or rank == 0:
    print ("set up in rank 0")
    pdb_lines = open("./1m2a.pdb","r").read()
    from LS49.sim.util_fmodel import gen_fmodel

    GF = gen_fmodel(resolution=10.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.7)
    A = GF.get_amplitudes()

    from LS49.spectra.generate_spectra import spectra_simulation
    SS = spectra_simulation()

    transmitted_info = dict(spectra = SS,
                            amplitudes = A)
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

  #for item,key in get_items(key=3271):
  for item,key in get_items(rank):
    result = dict(image=key,millers=[],orig_idx=[],pick=pickle_glob%key,spectrumx=[],obs=[],model=[],cc=[],
                  simtbx_millers=[],simtbx_intensity=[])
    #good indices: 3271, 2301: continue
    d = item["d"]
    nitem += 1

    print ("key %d in rank %d"%(key,rank))

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
        print ("key %d in rank %d hkl %s"%(key,rank,hkl[x]))
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
        values = sb.data-sb.background # ADU above background
        # intensity_lookup consists of the "observed" data from shoeboxes

        v0 = values.set_selected(values<=0, 0.)
        v1 = v0.set_selected(v0>255,255)
        v2 = (256.-v1)/256.
        np_v2 = np.ndarray(shape=(B[3]-B[2],B[1]-B[0],), dtype=np.float32, buffer=v2.as_numpy_array())

        # insert code here to estimate the partiality response
        PRD = get_partiality_response(key,hkl[x],spectra_simulation=transmitted_info["spectra"],ROI=ROI)
        pr_value=PRD["roi_pixels"]
        miller=PRD["miller"]
        intensity=PRD["intensity"]

        for c in range(nsb):
          intensity_lookup_1[(int(sb.coords()[c][1]),int(sb.coords()[c][0]))] = pr_value[c]
        assert len(intensity_lookup_1) == len(intensity_lookup)
        assert len(pr_value) == len(sb.data)
        # intensity_lookup_1 consists of partiality model data from posthoc simulator (partiality x Icalc)

        values_1 = pr_value # sb.data # ADU
        v0_1 = values_1.set_selected(values_1<=0, 0.)
        v1_1 = v0_1.set_selected(v0_1>255,255)
        v2_1 = (256.-v1_1)/256.
        np_v2_1 = np.ndarray(shape=(B[3]-B[2],B[1]-B[0],), dtype=np.float64, buffer=v2_1.as_numpy_array())

        ssplot=False
        if ssplot:
          ax=ax1=ax2=None
          fig = plt.figure(figsize=(8,6))
          ax = plt.subplot2grid ((2,2),(0,0))
          ax1 = plt.subplot2grid((2,2),(0,1))
          ax2 = plt.subplot2grid((2,2),(1,0),colspan=2)
          ax.imshow(np_v2, cmap=plt.cm.gray, interpolation="nearest")
          ax1.imshow(np_v2_1, cmap=plt.cm.gray, interpolation="nearest")
          plt.show()

        spectrumx,spectrumy,combined_model,CC =\
        plot_energy_scale_noplot(transmitted_info["spectra"],
                                 transmitted_info["amplitudes"].unit_cell().d(hkl[x]),
                                 abs_PA,origin,position0,B,intensity_lookup,intensity_lookup_1,key)
        result["simtbx_millers"].append(miller)
        result["simtbx_intensity"].append(intensity)
        result["millers"].append(asu[x])
        result["orig_idx"].append(hkl[x])
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
    import pickle
    print ("pickling key %d in rank %d"%(key,rank),result)
    pickle.dump(result,open("dataX%04d.pickle"%rank,"ab"))
