from __future__ import division, print_function
from six.moves import cPickle as pickle
from six.moves import range

with (open("debug26_range_intensities.pickle","rb")) as F:
    per_HKL_I = pickle.load(F)

with (open("debug26_intensities.pickle","rb")) as F:
    per_HKL_I_7122 = pickle.load(F)

from matplotlib import pyplot as plt


ten_hkl = [key for key in per_HKL_I][0:50]

for M in ten_hkl:
  print (M,per_HKL_I_7122[M])
  plt.plot(range(7090,7151),per_HKL_I[M])
  plt.plot([7122,], per_HKL_I_7122[M],"r.")

plt.show()
