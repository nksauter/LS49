from __future__ import division
import sys,os
logdir = sys.argv[1] or "."
allfiles = os.listdir(logdir)
basic_x=[]; renor_y=[]
for item in allfiles:
  ext = os.path.splitext(item)[1]
  if ext==".log":
    both = []
    with open(os.path.join(logdir,item),"r") as F:
      for line in F.readlines():
        tokens = line.strip().split()
        if "Basic" in tokens:
          both.append(float(tokens[7]))
        if "Renormalized" in tokens:
          both.append(float(tokens[4]))
        if len(both)==2:
          basic_x.append(both[0])
          renor_y.append(both[1])
          both = []

from matplotlib import pyplot as plt
plt.plot(basic_x, renor_y, "r.")
plt.plot([0.0,2.0],[0.0,2.0], "b--")
plt.title("Model vs. Experimental Bragg spot position")
plt.xlabel("DIALS rmsd (px)")
plt.ylabel("DiffBragg stage 1 rmsd (px)")
axes=plt.gca()
axes.set_aspect("equal")
plt.show()

