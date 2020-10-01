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
from matplotlib.ticker import FormatStrFormatter
from math import log,exp


fig, axes = plt.subplots(1,2,sharey=True,figsize=(12,6))


X = [60, 120, 240, 300]

js1 = [695, 399, 288, 216]
we1 = [657.4,355.4,239.5, 166.6]

js2 = [523, 295, 208, 204]
we2 = [489.0,249.6,147.3, 134.9]

js3 = [486, 284, 198, 206]
we3 = [448.8,228.7,126.7, 120]

js4 = [469, 272, 198, 203]
we4 = [427.7,220.0,124.5, 117.4]

js5 = [450, 267, 201, 209]
we5 = [410.3,215.1,121.6, 106.5]

js6 = [445, 268, 201, 202]
we6 = [401.0,215.0,113.4, 92.9]

js7 = [451, 262, 205, 219]
we7 = [398.6,204.5,114.1, 100.8]

colmap = ["#000000","#B30BFA","#0B1BFA","#0BBBFA","#0BFA33","#F2FA0B","#FA9B0B","#FA330B"]

ax = axes [0]

ax.set_title("""a) Strong scaling
on total execution time""")
ax.loglog(X, js1, color=colmap[1], marker="o")
ax.loglog(X, js2, color=colmap[2], marker="o")
ax.loglog(X, js3, color=colmap[3], marker="o")
ax.loglog(X, js4, color=colmap[4], marker="o")
ax.loglog(X, js5, color=colmap[5], marker="o")
ax.loglog(X, js7, color=colmap[7], marker="o")
ax.loglog(X, js6, color=colmap[6], marker="o")
loglogslope = (log(295.)-log(523.))/(log(120.)-log(60.))
extrap = exp( log(262) + loglogslope * (log(300)-log(120)) )
ax.loglog([120,300], [262, extrap],"--",color=colmap[7])
ax.set_xlim(50,330)
ax.set_ylim(80,800)
ax.xaxis.set_ticks(X, minor=True)
ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
ax.xaxis.set_ticks([], minor=False)
ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
ax.yaxis.set_ticks([100,200,300,400,500,600,700], minor=False)
ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
ax.set_xlabel("Number of nodes")
ax.set_ylabel("Total wall time (sec)")

ax = axes[1]
ax.set_title("""b) Strong scaling
on t > 0 GPU execution time""")
ax.loglog(X, we1, color=colmap[1], marker="o", label = "N=1")
ax.loglog(X, we2, color=colmap[2], marker="o", label = "N=2")
ax.loglog(X, we3, color=colmap[3], marker="o", label = "N=3")
ax.loglog(X, we4, color=colmap[4], marker="o", label = "N=4")
ax.loglog(X, we5, color=colmap[5], marker="o", label = "N=5")
ax.loglog(X, we6, color=colmap[6], marker="o", label = "N=6")
ax.loglog(X, we7, color=colmap[7], marker="o", label = "N=7")
ax.loglog(X, we6, color=colmap[6], marker="o")
ax.loglog(300,92.9, color=colmap[6], marker="o", markersize=15)
ax.set_xlim(50,330)
ax.set_ylim(80, 800)
ax.xaxis.set_ticks(X, minor=True)
ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
ax.xaxis.set_ticks([], minor=False)
ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
ax.yaxis.set_ticks([100,200,300,400,500,600,700], minor=False)
ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
#help (ax.legend)
ax.legend(loc="upper right",fontsize=12)
ax.set_xlabel("Number of nodes")

plt.show()
