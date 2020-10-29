from __future__ import division
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'xtick.major.width':   2})
matplotlib.rcParams.update({'ytick.major.width':   2})
matplotlib.rcParams.update({'xtick.minor.width':   2})
matplotlib.rcParams.update({'ytick.minor.width':   2})
matplotlib.rcParams.update({'xtick.major.size':  8})
matplotlib.rcParams.update({'ytick.major.size':   8})
matplotlib.rcParams.update({'xtick.minor.size':   8})
matplotlib.rcParams.update({'ytick.minor.size':   0})
matplotlib.rcParams.update({'axes.linewidth':   2})
matplotlib.rcParams.update({'axes.grid':       True})
matplotlib.rcParams.update({'axes.grid':       True})
matplotlib.rcParams.update({'axes.grid.which': "both"})
matplotlib.rcParams.update({'grid.alpha': 0.5})
from matplotlib.ticker import FormatStrFormatter

fig, axes = plt.subplots(1,1,sharey=True,figsize=(7,6))

X = [1, 2, 4, 7, 15, 30, 60, 119, 238, 476]

we1 = [128.5, 129.6, 130.7, 133.0, 132.7, 132.9, 132.7, 133.2, 152.5, 153.1]

we2 = [166, 169, 170, 167, 172, 175, 187, 182, 206, 211]

we1_1device_only = [684.1, 684.4, 684.3, 691.6, 689.9, 688.8, 694.7, 694.0, 794.2, 803.7]

colmap = ["#000000","#B30BFA","#0B1BFA","#0BBBFA","#0BFA33","#F2FA0B","#FA9B0B","#FA330B"]

ax = axes
one_device_only = True
if one_device_only:
  ax.loglog(X, we1_1device_only, color=colmap[4], marker="o", label = "1 GPU")
ax.loglog(X, we2, color=colmap[2], marker="o", label = "setup + 6 GPUs")
ax.loglog(X, we1, color=colmap[7], marker="o", label = "6 GPUs")


#ax.loglog(300,91.3, color=colmap[7], marker="o", markersize=15)
ax.set_xlim(0.75,575)
ax.set_ylim(80, 900)

Xax = [1,2,4,8,16,32,64,128,256,512]
ax.xaxis.set_ticks(Xax, minor=True)
ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
ax.xaxis.set_ticks([], minor=False)
ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
ax.yaxis.set_ticks([100,200,300,400,500,600,700,800], minor=False)
ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
#help (ax.legend)
ax.legend(loc="center right",fontsize=12)
ax.set_xlabel("Number of nodes")
ax.set_ylabel("Total wall time (sec)")

plt.show()
