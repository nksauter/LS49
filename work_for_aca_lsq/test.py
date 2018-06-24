from __future__ import division, print_function
from six.moves import range
import dxtbx

#file1 = "step5_000000.img.gz"
#file2 = "../../../LS49/step5_000000.img"

file1 = "../../../LS49/step5Bpoly_000000.img.gz"
file2 = "../../../LS49/step5_MPIbatch_000000.img.gz"

file1 = "step5poly_000000.img.gz"
file2 = "../test1/step5poly_000000.img.gz"

L1 = dxtbx.load(file1)
L2 = dxtbx.load(file2)

d1 = L1.get_raw_data()
d2 = L2.get_raw_data()

for i in range(len(d1)):
 if d1[i]-d2[i] != 0:
  print("%6d %6d %6d"%(i, d1[i], d2[i]), end=' ')
  if d1[i]-d2[i] != 0:
    print(d1[i]-d2[i])
  else:
    print()
