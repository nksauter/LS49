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
from math import log,exp


fig, axes = plt.subplots(1,2,sharey=True,figsize=(13,6))


X = [60, 120, 240, 300]

Tera_ops = [
 0.01718,0.03221,0.06441,0.13312,0.2684,0.5368,0.8052,1.0735,1.6103,2.1471,3.2206]
Mullen_argchk_addspots_wall = [
 1.25,1.89,3.04,5.81,11.28,22.04,32.94,43.63,65.46,91.04,138.42]
Debranched_addspots_wall = [
 0.92,1.33,2.27,3.75,7.33,14.39,21.55,28.53,42.74,67.57,87.01]


colmap = ["#000000","#B30BFA","#0B1BFA","#0BBBFA","#0BFA33","#F2FA0B","#FA9B0B","#FA330B"]

ax = axes[0]

ax.set_title("""a) Jungfrau, adse13-187""")
loglogslope = 1.
extrap = exp( log(53.31) - loglogslope * (log(2.88)-log(0.0225)) )
ax.loglog([1.0735,], [43.63,], color=colmap[2], marker="o", markersize=10) # emphasis canonical problem size
ax.loglog([0.0225,2.88], [extrap,53.31],"--",color=colmap[7], label="ideal scaling")

ax.loglog(Tera_ops, Mullen_argchk_addspots_wall, color=colmap[2], marker="o",
          label="Mullen kernel")
ax.loglog(Tera_ops, Debranched_addspots_wall, color=colmap[4], marker="o",
          label="Debranched kernel")

ax.set_xlim(0.014,3.6)
#ax.set_ylim(80,800)
ax.xaxis.set_ticks([0.1,1], minor=False)
ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
ax.yaxis.set_ticks([1,10,100], minor=False)
ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
ax.set_xlabel("Number of arg evaluations (Tera-ops)")
ax.set_ylabel("Bragg spots wall time (sec)")

ax.legend(loc="upper left",fontsize=12)
Tera_ops = [0.0225,0.045,0.09,0.18,0.36,0.72,1.08,1.44,2.16,2.88]
Mullen_argchk_addspots_wall = [1.52,2.16,3.55,6.4,12.14,20.23,30.14,48.85,76.05,87.23]
Debranched_addspots_wall = [0.88,1.25,2.02,3.64,6.91,13.45,20.01,26.56,40.11,53.31]

ax = axes[1]

ax.set_title("""b) LS49, adse13-196""")
loglogslope = 1.
extrap = exp( log(53.31) - loglogslope * (log(2.88)-log(0.0225)) )
ax.loglog([0.0225,], [1.52,], color=colmap[2], marker="o", markersize=10) # emphasis canonical problem size
ax.loglog([0.0225,2.88], [extrap,53.31],"--",color=colmap[7], label="ideal scaling")

ax.loglog(Tera_ops, Mullen_argchk_addspots_wall, color=colmap[2], marker="o",
          label="Mullen kernel")
ax.loglog(Tera_ops, Debranched_addspots_wall, color=colmap[4], marker="o",
          label="Debranched kernel")

ax.set_xlim(0.014,3.6)
#ax.set_ylim(80,800)
ax.xaxis.set_ticks([0.1,1], minor=False)
ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
ax.yaxis.set_ticks([1,10,100], minor=False)
ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
ax.set_xlabel("Number of arg evaluations (Tera-ops)")

ax.legend(loc="upper left",fontsize=12)

plt.show()
