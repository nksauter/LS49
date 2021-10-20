from __future__ import division
comb1  = "/global/cscratch1/sd/nksauter/adse13_187/bleededge/work/2310293/combined.refl"
comb2  = "/global/cscratch1/sd/nksauter/adse13_187/bleededge/work/residuals/combined.refl"
from dials.array_family import flex

def proc(comb):
  refl = flex.reflection_table.from_file(comb)
  diff_vec3 = refl["xyzcal.px"] - refl["xyzobs.px.value"]
  sqdiff = diff_vec3.dot(diff_vec3)
  beam_vec3 = flex.vec3_double(len(refl), [0.,0.,-1.])
  dotbeam = (refl["s1"].each_normalize()).dot(beam_vec3.each_normalize())
  angle = flex.acos(dotbeam)
  refl["angle"]=angle
  refl["sqdiff"]=sqdiff
  print(diff_vec3,sqdiff)
  print("%d refls"%len(refl))
  order = flex.sort_permutation(angle)
  refl["order"]=order
  return refl

def x_y(r):
  sortr = r["sqdiff"].select(r["order"]) # square diffs, sorted
  xvals = [500 + 1000*N for N in range(0, len(r)//1000)]
  yvals = []
  for N in range(0, len(r)//1000):
    interval = sortr[1000*N:1000*(N+1)]
    mean = flex.mean(interval)
    import math
    rms = math.sqrt(mean); yvals.append(rms)
  assert len(xvals)==len(yvals)
  return xvals,yvals

r1 = proc(comb1)
r2 = proc(comb2)
x1,y1 = x_y(r1)
x2,y2 = x_y(r2)
from matplotlib import pyplot as plt
#plt.plot([1,2,3],[4,5,6])
#plt.plot(r1["angle"],r1["sqdiff"],"r.")
#plt.plot(r2["angle"],r2["sqdiff"],"b.")
plt.plot(x1,y1,"r.")
plt.plot(x2,y2,"b.")
plt.plot(x1,y1,"r-")
plt.plot(x2,y2,"b-")
plt.show()
