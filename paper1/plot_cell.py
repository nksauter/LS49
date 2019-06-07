from __future__ import division, print_function
from scitbx.array_family import flex
# dials.show idx-step5_MPIbatch_00[0-5]*.img_integrated_experiments.json|grep "Unit cell" > cellrestrain.5000
# dials.show ../LS49_integ_step5cori/idx-step5_MPIbatch_00[0-5]*.img_integrated_experiments.json|grep "Unit cell" > cellnorestrain.5000
params_rest = [flex.double(), flex.double(), flex.double(), flex.double()]
items = ["a","b","c","beta"]

with open("cellrestrain.10000") as P:
  for ili,line in enumerate(P):
    tokens = line.strip().split()
    a = float (tokens[2].split("(")[1])
    b = float (tokens[3].split("(")[0])
    c = float (tokens[4].split("(")[0])
    beta = float(tokens[6].split("(")[0])
    for idx,item in enumerate([a,b,c,beta]):  params_rest[idx].append(item)
for idx in range(4):
  stats = flex.mean_and_variance(params_rest[idx])
  print(items[idx],"restrained mean",stats.mean(),
        "sigma",stats.unweighted_sample_standard_deviation(),
        "on",len(params_rest[idx])
       )

paramsnorest = [flex.double(), flex.double(), flex.double(), flex.double()]
with open("cellnorestrain.10000") as P:
  for line in P:
    tokens = line.strip().split()
    try:
      a = float (tokens[2].split("(")[1])
      b = float (tokens[3].split("(")[0])
      c = float (tokens[4].split("(")[0])
      beta = float(tokens[6].split("(")[0])
      for idx,item in enumerate([a,b,c,beta]):  paramsnorest[idx].append(item)
    except Exception: pass
for idx in range(4):
  stats = flex.mean_and_variance(paramsnorest[idx])
  print(items[idx],"unrestrained mean",stats.mean(),
        "sigma",stats.unweighted_sample_standard_deviation(),
        "on",len(paramsnorest[idx])
       )

from matplotlib import pyplot as plt
fig, axes = plt.subplots( nrows=1, ncols=4, figsize=(10,6), sharey=True )

for idx, leg in enumerate(["cell a", "cell b", "cell c", "beta angle"]):
  axes[idx].set_title(leg)
n,bins,patches = axes[3].hist(paramsnorest[3],40, range=(110.115,110.515), normed=0, facecolor="blue", alpha=0.75)
n,bins,patches = axes[3].hist(params_rest[3],40, range=(110.115,110.515), normed=0, facecolor="orange", alpha=0.75)

n,bins,patches = axes[0].hist(paramsnorest[0],40, range=(67.0,67.4), normed=0, facecolor="blue", alpha=0.75)
n,bins,patches = axes[0].hist(params_rest[0],40, range=(67.0,67.4), normed=0, facecolor="orange", alpha=0.75)

n,bins,patches = axes[1].hist(paramsnorest[1],40, range=(59.6,60.0), normed=0, facecolor="blue", alpha=0.75)
n,bins,patches = axes[1].hist(params_rest[1],40, range=(59.6,60.0), normed=0, facecolor="orange", alpha=0.75)

n,bins,patches = axes[2].hist(paramsnorest[2],40, range=(47.05,47.35), normed=0, facecolor="blue", alpha=0.75)
n,bins,patches = axes[2].hist(params_rest[2],40, range=(47.05,47.35), normed=0, facecolor="orange", alpha=0.75)
plt.show()
