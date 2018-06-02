from __future__ import division, absolute_import
from six.moves import range
from post5_ang_misset import parse_postrefine
from scitbx.matrix import col
from dials.algorithms.shoebox import MaskCode
from scitbx.array_family import flex
json_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/idx-step5_MPIbatch_0%05d.img_integrated_experiments.json"
image_glob = "/net/dials/raid1/sauter/LS49/step5_MPIbatch_0%05d.img.gz"
pickle_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/idx-step5_MPIbatch_0%05d.img_integrated.pickle"

pdb_lines = open("/net/dials/raid1/sauter/LS49/1m2a.pdb","r").read()
from LS49.sim.util_fmodel import gen_fmodel
GF = gen_fmodel(resolution=10.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.7)
A = GF.get_amplitudes()

def get_items():
  for key in postreffed:
    print key
    from dxtbx.model.experiment_list import ExperimentListFactory
    E = ExperimentListFactory.from_json_file(json_glob%key,check_format=False)[0]
    C = E.crystal
    from six.moves import cPickle as pickle
    T = pickle.load(open(pickle_glob%key,"rb"))
    #from IPython import embed; embed()
    resolutions = T["d"]
    millers = T["miller_index"]
    nitem = len(resolutions)
    yield T


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
  for item in get_items():
    d = item["d"]
    nitem += 1
    #if nitem>29525:break
    nall_spots += len(item)
    iselect = ((d < 2.5) & (d > 2.1))
    nres_range += len(d.select(iselect))

    # geometric selection:  between position angles 150 and 210 degrees.
    hkl = item["miller_index"].select(iselect)
    cust_copy = A.customized_copy(indices=hkl,data=flex.double(len(hkl)),sigmas=flex.double(len(hkl)))
    asu = cust_copy.map_to_asu().indices()
    xyz = item["xyzobs.px.value"].select(iselect)
    shoe = item["shoebox"].select(iselect)
    for x in range(len(hkl)):
      slow = xyz[x][1]
      fast = xyz[x][0]
      positionX = col((slow,fast))-origin
      position_angle = positionX.angle(position0,deg=True)
      if position_angle > 150.:
        npos_angle += 1
        millerd[asu[x]]=millerd.get(asu[x],0)+1
        print "%20s,asu %20s"%(str(hkl[x]),str(asu[x])),"slow=%5.0f  fast=%5.0f"%(xyz[x][1],xyz[x][0]),"PA %6.1f"%position_angle
        #from IPython import embed; embed()
        sb = shoe[x]
        #print "shoebox has %d pixels"%(sb.mask.size())
        nsb = sb.mask.size()
        for c in range(nsb):
          #print c, sb.mask[c], sb.mask[c]&MaskCode.Valid == MaskCode.Valid,\
          #         sb.mask[c]&MaskCode.Foreground == MaskCode.Foreground,\
          #         sb.mask[c]&MaskCode.Background == MaskCode.Background,\
          #         sb.mask[c]&MaskCode.BackgroundUsed   == MaskCode.BackgroundUsed
          if sb.mask[c]&MaskCode.Valid == MaskCode.Valid and sb.mask[c]&MaskCode.Foreground == MaskCode.Foreground:
            #print "VF",c, sb.coords()[c], sb.background[c], sb.data[c], sb.data[c]-sb.background[c]
            nVF += 1


  print "Number of images %d; of all spots %d; of in-resolution spots %d; in position %d"%(
    nitem, nall_spots, nres_range, npos_angle)
  print "Valid foreground pixels: %d. Number of Miller indices: %d"%(nVF, len(millerd))
  print "Average",  npos_angle/nitem,"spots/image"
  print "Average",  npos_angle/len(millerd), "observations/Miller index"
  print "Average",  nVF/npos_angle," valid foreground pixels /spot"
  nfreq = 0
  print "Analyze Miller indices observed more than 30 times"
  for key in millerd:
    if millerd[key]>30:
      print key, millerd[key]
      nfreq+=1
  print "Total of %d observed > 30 times"%nfreq
