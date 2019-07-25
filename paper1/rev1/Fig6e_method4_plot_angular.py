from __future__ import division, absolute_import, print_function
import sys,math
trial = "slurm22693509.out" # sys.argv[1] # Like "005"
histocolor = "red" # sys.argv[2] # blue for method 1, orange for method 2
from scitbx.array_family import flex
if True:
  angular_offset = flex.double()


  for line in open(trial,"r"):
    if line.find("LLG")>=0:
      value = line.split()[19]
      angular_offset.append(float(value))

  # ersatz selection to get rid of 18 misindexed outliers in Method 1
  angular_offset = angular_offset.select( angular_offset < 37. )

  order = flex.sort_permutation(angular_offset)
  sorted_missettings = angular_offset.select(order)
  print(list(sorted_missettings))
  print ("Count > 0.04 degrees",(angular_offset>0.04).count(True),"of total",len(angular_offset))


  median = angular_offset[order[len(order)//2]]
  rmsd = math.sqrt(flex.mean(angular_offset*angular_offset))
  print(trial,"%6d measurements; rmsd %9.5f"%(len(angular_offset),
    rmsd
  ),
    "Median is %9.5f"%(median))
  print("Max is %9.5f"%(flex.max(angular_offset)))

  from matplotlib import pyplot as plt
  zoom="hi"
  #nbins = dict(hi=300,lo=100)[zoom]
  # to ensure bins always represent 0.001, multiply mx by 1000
  nbins = int(1000.*flex.max(angular_offset))
  n,bins,patches = plt.hist(angular_offset,
    nbins, normed=0, facecolor=histocolor, alpha=0.75) # orange or blue

  plt.xlabel("Orientational offset from truth (degrees)")
  plt.title("Histogram of cctbx.xfel misorientation %s"%trial)
  maxangle = 0.15
  PLOT_Height = 11800
  plt.axis([-0.001,maxangle,0,PLOT_Height])

  plt.annotate(
        "median %6.4f"%median,
        xy=(median, 0.2*PLOT_Height), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
  plt.annotate(
        "rmsd %6.4f"%rmsd,
        xy=(rmsd, 0.325*PLOT_Height), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

  mean = flex.mean(angular_offset)
  plt.annotate(
        "mean %6.4f"%mean,
        xy=(mean, 0.45*PLOT_Height), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

  plt.plot([median,rmsd,mean],[0.2*PLOT_Height,0.325*PLOT_Height,0.45*PLOT_Height],"b|")
  plt.show()
