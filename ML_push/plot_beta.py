from __future__ import division, print_function
from scitbx.array_family import flex
# cat ~/ML_push/merge5_redo2.digest | grep -A1 integrated | grep Unit > cells.10000
betas = flex.double()
with open("cells.10000") as P:
  for line in P:
    tokens = line.strip().split()
    beta = float(tokens[6].replace(",",""))
    betas.append(beta)
stats = flex.mean_and_variance(betas)
print("mean",stats.mean(),"sigma",stats.unweighted_sample_standard_deviation())
from matplotlib import pyplot as plt
n,bins,patches = plt.hist(betas,40, range=(110.115,110.515), normed=0, facecolor="blue", alpha=0.75)
plt.show()
