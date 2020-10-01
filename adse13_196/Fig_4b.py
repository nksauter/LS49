from __future__ import division
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'xtick.major.width':   2})
matplotlib.rcParams.update({'ytick.major.width':   2})
matplotlib.rcParams.update({'axes.linewidth':   2})

X = [1,2,3,4,5,6,7]

wtime = [2.34,3.45,4.76,5.96,7.22,8.43,9.62]

unit_time = [2.34, 1.73, 1.59, 1.49, 1.44, 1.41, 1.37]

speedup = [ 2.34 * 7 / (wtime[i] * (7./X[i])) for i in range(7)]


import numpy as np
slope,y_intercept = np.polyfit( X, wtime, 1)
print (slope, y_intercept, "params")

for x,y in zip(X,wtime):
  print (x, y, slope*x+y_intercept, speedup[x-1])

plt.plot(X, wtime, "bo")
plt.plot(X, wtime, "b-")
#plt.plot(X, unit_time, "go")
#plt.plot(X, unit_time, "g-")
#plt.plot(X, speedup,"r-")

plt.xlim(0,8)
plt.ylim(0,10)

plt.show()
