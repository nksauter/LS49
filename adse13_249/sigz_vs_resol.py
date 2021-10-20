from __future__ import division
comb1  = "combined.refl"
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
expts = ExperimentListFactory.from_json_file("combined.expt",
                                              check_format=False)

def proc(comb,expts):
  refl = flex.reflection_table.from_file(comb)
  beam_vec3 = flex.vec3_double(len(refl), [0.,0.,-1.])
  dotbeam = (refl["s1"].each_normalize()).dot(beam_vec3.each_normalize())
  angle = flex.acos(dotbeam)
  refl["angle"]=angle
  refl["d"] = refl.compute_d_single(expts[0]) # sort by this rather than angle, show Miller
  refl["dstar"] = 1./refl["d"]
  order = flex.sort_permutation(refl["dstar"])
  refl["order"]=order
  print(list(refl["order"]))
  #print ( list(refl["dstar"].select(refl["order"]))  )
  return refl

class binner:
  def __init__(self, refl, n):
    self.refl = refl
    self.n = n
    dstar_max = flex.max(refl["dstar"])
    self.recip_volume_series=flex.double([(dstar_max**3)*((x+1.0)/self.n) for x in range(self.n)])
    self.recip_volume_labels=flex.double([(dstar_max**3)*((x+0.5)/self.n) for x in range(self.n)])
    self.d_series = flex.pow( self.recip_volume_series, -1./3. )
    self.d_labels = flex.pow( self.recip_volume_labels, -1./3. )
    self.refl["dstar3"] = flex.pow(self.refl["dstar"], 3)
    dstar3_max = flex.max(refl["dstar3"])
    print(list(self.d_series))
    print(list(self.d_labels))
    self.bin_selections=[]
    for x in range(self.n):
      self.bin_selections.append(
        (self.refl["dstar3"] > (self.recip_volume_series[x-1] if x>0 else 0.)
        ).__and__(
        self.refl["dstar3"] <= (self.recip_volume_series[x]))
        )
    self.bin = flex.int( [int(x) for x in  ((self.refl["dstar3"]/dstar3_max)*self.n)-1.0] )
    print('len', len(self.bin))
    print(list(self.bin))
    self.sort_bin = self.bin.select(self.refl["order"])
    self.sort_sgz = self.refl["sigma_z"].select(self.refl["order"])
    self.bincts = [sel.count(True) for sel in self.bin_selections]
  def ranges(self):
    ptr = 0
    for x in range(self.n):
      xvals = range(ptr,ptr+self.bincts[x])
      yvals = self.sort_sgz[ptr:ptr+self.bincts[x]]
      ptr += self.bincts[x]
      yield xvals,yvals

def conv(r):
  return range(len(r)), r["sigma_z"].select(r["order"])

def x_y(r):
  sortr = r["sigma_z"].select(r["order"])
  print(["%.1f"%x for x in list(sortr)])
  xvals = [5 + 10*N for N in range(0, len(r)//10)]
  yvals = []
  for N in range(0, len(r)//10):
    interval = sortr[10*N:10*(N+1)]
    mean = flex.mean(interval)
    import math
    rms = math.sqrt(mean); yvals.append(rms)
  assert len(xvals)==len(yvals)
  return xvals,yvals

r1 = proc(comb1,expts)
B = binner(r1,n=10)
from matplotlib import pyplot as plt

def pltall():
    x1,y1 = conv(r1)#x_y(r1)
    plt.plot(x1,y1,"r-")
    plt.show()
#pltall()

def pltbins():
  for x,y in B.ranges():
    plt.plot(x,y,".")
  plt.show()

def pltstats():
  print("\nQuartile distributions of σ<Z> for %d reflections on %d lattices"%(
    len(r1), len(expts)))
  print("Bin   Resolution    Nrefl median(σ<Z>)      Q1-Q3")
  fig, axes = plt.subplots(1,1)
  ax = axes # only one
  kvals = []
  ivals = []
  for i,(x,y) in enumerate(B.ranges()):
    if len(y)<4: continue
    S = sorted(y)
    median = S[len(y)//2]
    Q1 = S[len(y)//4]
    Q3 = S[3*len(y)//4]
    plt.plot(i, median,"b.")
    plt.plot([i,i],[Q1,Q3],"b-")
    kvals.append(median)
    ivals.append(i)
    print(" %2d  %5s -%5.2f  %6d    %.2f          %.2f - %.2f"%(
     i, "%5.2f"%(B.d_series[i-1]) if i>0 else "   ∞ ",
     B.d_series[i], len(S), median, Q1, Q3
    ))
  plt.plot(ivals, kvals,"r-")
  ax.set_title("Quartile distributions of σ<Z> for\n%d reflections on %d lattices"%(
    len(r1), len(expts)))
  ax.set_ylabel("Shoebox σ<Z> (std. deviations)")
  ax.set_xlabel("Resolution (Å)")
  ax.xaxis.set_ticks(range(B.n), minor=False)
  ax.xaxis.set_ticklabels(["%.2f"%x for x in B.d_labels])
  plt.show()
pltstats()

exit()
for selection in B.bin_selections:
  print(list(selection))
  d = r1["d"].select(selection)
  order = r1["order"].iselect(selection)
  sigmaz = r1["sigma_z"].select(selection)
  print(len(order), len(sigmaz))
  plt.plot(order, sigmaz, ".")
  break
plt.show()

#print(["%.1f"%x for x in list(r1["sigma_z"])])
exit()
x1,y1 = conv(r1)#x_y(r1)
plt.plot(x1,y1,"r.")
plt.plot(x1,y1,"r-")
plt.show()
