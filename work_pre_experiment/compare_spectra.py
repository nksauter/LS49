from post5_ang_misset import parse_postrefine
from scitbx.matrix import col
from dials.algorithms.shoebox import MaskCode
from scitbx.array_family import flex
from matplotlib.ticker import FuncFormatter
import numpy as np
import math

json_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/idx-step5_MPIbatch_0%05d.img_integrated_experiments.json"
image_glob = "/net/dials/raid1/sauter/LS49/step5_MPIbatch_0%05d.img.gz"
pickle_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/idx-step5_MPIbatch_0%05d.img_integrated.pickle"

pdb_lines = open("/net/dials/raid1/sauter/LS49/1m2a.pdb","r").read()
from LS49.sim.util_fmodel import gen_fmodel
GF = gen_fmodel(resolution=10.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.7)
A = GF.get_amplitudes()

#specialize this file to look at one particular index
one_index = (14,-22,-2)
distance_mm = 141.7
pixel_sz_mm = 0.11
mos_rotation_deg = 0.05

from LS49.spectra.generate_spectra import spectra_simulation
SS = spectra_simulation()

def get_items():
  for key in postreffed:
    print "image=",key
    from dxtbx.model.experiment_list import ExperimentListFactory
    E = ExperimentListFactory.from_json_file(json_glob%key,check_format=False)[0]
    C = E.crystal
    import cPickle as pickle
    T = pickle.load(open(pickle_glob%key,"rb"))
    #from IPython import embed; embed()
    resolutions = T["d"]
    millers = T["miller_index"]
    nitem = len(resolutions)
    yield T,key

def plot_energy_scale(ax,ax2,abs_PA,origin,position0,B,intensity_lookup,key):
  d_Ang = A.unit_cell().d(one_index)
  print d_Ang
  unit_pos0 = position0.normalize()
  spectrumx = []
  spectrumy = []
  for eV in range(7090,7151):
    spectrumx.append(eV)
    specy = 0.
    lambda_Ang = 12398.425 / eV
    two_theta = 2. * math.asin( lambda_Ang / (2.*d_Ang))
    radius_mm = distance_mm * math.tan(two_theta)
    radius_px = radius_mm / pixel_sz_mm
    contour_x = []; contour_y = []
    for rot in range(-2,3):
      PA = abs_PA + rot*mos_rotation_deg*math.pi/180.
      clock = unit_pos0.rotate_2d(-PA, deg=False)
      position1 = origin + radius_px*clock
      contour_x.append(position1[1]-0.5-B[0])
      contour_y.append(position1[0]-0.5-B[2])
    
    ax.plot(contour_x,contour_y, "r-",linewidth=0.3)

    contour_x = []; contour_y = []
    for rot in range(-8,9):
      PA = abs_PA + 0.25*rot*mos_rotation_deg*math.pi/180.
      clock = unit_pos0.rotate_2d(-PA, deg=False)
      position1 = origin + radius_px*clock
      int_coords = (int(position1[0]),int(position1[1]))
     
      specy += intensity_lookup.get(int_coords,0)
      if int(position1[1])!=367 : continue
      if int(position1[0])!=1075: continue
      contour_x.append(position1[1]-0.5-B[0])
      contour_y.append(position1[0]-0.5-B[2])
    
    ax.plot(contour_x,contour_y, "ro",markersize=1.0)
    spectrumy.append(specy)
  ax2.plot(spectrumx, spectrumy,"r.")
  iterator = SS.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)
  wavlen, flux, wavelength_A = iterator.next() # list of lambdas, list of fluxes, average wavelength
  ratio = flex.max(flux)/max(spectrumy)
  ax2.plot(12398.425/wavlen,(flux/ratio)+100,"b-") 

if __name__=="__main__":
  origin = col((1500,1500))
  position0 = col((1500,3000))-origin
  postreffed = parse_postrefine()
  print "# postrefined images",len(postreffed)
  nitem = 0
  nall_spots = 0
  nres_range = 0
  npos_angle = 0
  nVF = 0
  millerd = {}
  for item,key in get_items():
    d = item["d"]
    nitem += 1
    
    nall_spots += len(item)
    iselect = ((d < 2.5) & (d > 2.1))
    nres_range += len(d.select(iselect))

    # geometric selection:  between position angles 150 and 210 degrees.
    hkl = item["miller_index"].select(iselect)
    cust_copy = A.customized_copy(indices=hkl,data=flex.double(len(hkl)),sigmas=flex.double(len(hkl)))
    asu = cust_copy.map_to_asu().indices()
    xyz = item["xyzobs.px.value"].select(iselect)
    calcpx = item["xyzcal.px"].select(iselect)
    shoe = item["shoebox"].select(iselect)
    intensity_lookup ={}
    for x in xrange(len(hkl)):
      if not asu[x]==one_index:continue
      slow = xyz[x][1]
      fast = xyz[x][0]
      positionX = col((slow,fast))-origin
      position_angle = positionX.angle(position0,deg=True)
      if position_angle > 150.:
        npos_angle += 1
        millerd[asu[x]]=millerd.get(asu[x],0)+1
        print "%20s,asu %20s"%(str(hkl[x]),str(asu[x])),
        print "slow=%5.0f  fast=%5.0f"%(xyz[x][1],xyz[x][0]),"PA %6.1f"%position_angle
        sb = shoe[x]
        print "shoebox has %d pixels"%(sb.mask.size())  
        print "first coordinate is ",sb.coords()[0]
        nsb = sb.mask.size()
        for c in range(nsb):
          #print c, sb.mask[c], sb.mask[c]&MaskCode.Valid == MaskCode.Valid,\
          #         sb.mask[c]&MaskCode.Foreground == MaskCode.Foreground,\
          #         sb.mask[c]&MaskCode.Background == MaskCode.Background,\
          #         sb.mask[c]&MaskCode.BackgroundUsed   == MaskCode.BackgroundUsed
          if sb.mask[c]&MaskCode.Valid == MaskCode.Valid and sb.mask[c]&MaskCode.Foreground == MaskCode.Foreground:
            print "VF",c, sb.coords()[c], sb.background[c], sb.data[c], sb.data[c]-sb.background[c]
            nVF += 1
          intensity_lookup[(int(sb.coords()[c][1]),int(sb.coords()[c][0]))] = sb.data[c]-sb.background[c]
        print "bbox",sb.bbox
        print "spotprediction", calcpx[x]
        spotprediction = calcpx[x] # DIALS coords (fast,slow)
        spotvec = col((spotprediction[0],spotprediction[1]))-origin # panel coords (center=origin)
        print "spotvec",spotvec.elems
        abs_PA = math.atan2(spotvec[1],spotvec[0]) # clockwise plotted on image_viewer (because vertical axis upside down)
        print "abs position angle",abs_PA*180./math.pi
        B = sb.bbox
        values = sb.data-sb.background # ADU above background
        #from IPython import embed; embed()
        #for c in range(nsb):
        #  if (sb.mask[c]&MaskCode.Valid == MaskCode.Valid and sb.mask[c]&MaskCode.Foreground == MaskCode.Foreground):
        #    values[c]=0.
        v0 = values.set_selected(values<=0, 0.)
        v1 = v0.set_selected(v0>255,255)
        v2 = (256.-v1)/256.
        np_v2 = np.ndarray(shape=(B[3]-B[2],B[1]-B[0],), dtype=np.float32, buffer=v2.as_numpy_array())
        import matplotlib.pyplot as plt
        fig, (ax,ax2) = plt.subplots(2,1)
        ax.imshow(np_v2, cmap=plt.cm.gray, interpolation="nearest")

        def slow_scale(x, pos):
          return "%.1f"%(x + 0.5 + B[2])
        formatters = FuncFormatter(slow_scale)
        ax.yaxis.set_major_formatter(formatters)

        def fast_scale(x, pos):
          return "%.1f"%(x + 0.5 + B[0])
        formatterf = FuncFormatter(fast_scale)
        ax.xaxis.set_major_formatter(formatterf)
        
        ax.plot([spotprediction[0] - 0.5 - B[0]], [spotprediction[1] - 0.5 - B[2]], "y.")

        plot_energy_scale(ax,ax2,abs_PA,origin,position0,B,intensity_lookup,key)

        plt.show()

    
  print "Number of images %d; of all spots %d; of in-resolution spots %d; in position %d"%(
    nitem, nall_spots, nres_range, npos_angle)
  print "Valid foreground pixels: %d. Number of Miller indices: %d"%(nVF, len(millerd.keys()))
  print "Average",  npos_angle/nitem,"spots/image"
  print "Average",  npos_angle/len(millerd.keys()), "observations/Miller index"
  print "Average",  nVF/npos_angle," valid foreground pixels /spot"
  nfreq = 0
  print "Analyze Miller indices observed more than 30 times"
  for key in millerd:
    if millerd[key]>30:
      print key, millerd[key]
      nfreq+=1
  print "Total of %d observed > 30 times"%nfreq





    
