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

js1 = [729, 385, 234, 192]
we1 = [686.5,342.4,188.2,148.4]

js2 = [545, 291, 178, 142]
we2 = [502.0,251.7,143.1,102.8]

js3 = [504, 275, 173, 140]
we3 = [460.0,232.1,133.6,100.4]

js4 = [541, 269, 165, 144]
we4 = [502.2,226.6,124.2, 93.6]

js5 = [518, 262, 163, 140]
we5 = [479.0,220.5,124.0, 99.9]

js6 = [457, 265, 153, 142]
we6 = [414.6,221.7,107.2, 95.9]

js7 = [458, 260, 170, 141]
we7 = [415.1,211.7,122.9, 91.3]

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
loglogslope = (log(260.)-log(458.))/(log(120.)-log(60.))
extrap = exp( log(260.) + loglogslope * (log(300)-log(120)) )
ax.loglog([120,300], [260., extrap],"--",color=colmap[7])
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
loglogslope = (log(211.7)-log(415.1))/(log(120.)-log(60.))
extrap = exp( log(211.7) + loglogslope * (log(300)-log(120)) )
ax.loglog([120,300], [211.7, extrap],"--",color=colmap[7])


ax.loglog(300,91.3, color=colmap[7], marker="o", markersize=15)
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
