from __future__ import division
import sys,os
basic_x=[]; renor_y=[]; idx=[]; wt=[]
def gen1(logdir = sys.argv[1] or "."):
  allfiles = os.listdir(logdir)
  for item in allfiles:
    ext = os.path.splitext(item)[1]
    if ext==".log":
      both = []
      with open(os.path.join(logdir,item),"r") as F:
        for line in F.readlines():
          tokens = line.strip().split()
          if "Basic" in tokens:
            both = [float(tokens[8])]
          if "Renormalized" in tokens:
            both.append(float(tokens[6]))
            both.append(float(tokens[9])) # the number of shoeboxes on image
          if len(both)==3:
            basic_x.append(both[0])
            renor_y.append(both[1])
            wt.append(both[2]/10.)
            both = []
def gen2():
  for logdir in ["2281386","2281387","2281394","2281395","2281396", # 1 node, 20 min, 80 images
                 "2281498","2281499","2281523","2281524",           # 5 node, 20 min, 400
                 "2281545","2281546","2281547","2281548","2281549", # "
                 "2281557","2281558","2281560","2281561","2281562", # "
                 "2281570",
                 # 5 node, 80 min, 1600 (actually 1534) images.  All with --exclusive flag
                ]:
    gen1(logdir=logdir)
gen1()
#gen2()
print("analyzed",len(basic_x),"lattices")
from matplotlib import pyplot as plt
#plt.plot(basic_x,renor_y, "r.")
plt.scatter(basic_x,renor_y, s=wt, c="r", marker=".")
plt.plot([0.0,2.0],[0.0,2.0], "b--")
plt.title("Model vs. Experimental Bragg spot position")
plt.xlabel("DIALS rmsd (px)")
plt.ylabel("DiffBragg stage 1 rmsd (px)")
axes=plt.gca()
axes.set_aspect("equal")
plt.show()

