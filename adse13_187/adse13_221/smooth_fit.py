from __future__ import division
import numpy as np

def replacement_pixels(self, ipanel, islow_limits, ifast_limits):
  """main idea: LUNUS marks invalid model pixels as -1.0.  If so, try to find an equal number of valid
  pixels from around the edges of the box, prior to fitting the 2nd-order Taylor function.
  """
  A = all_lunus = self.lunus_filtered_data[ipanel, islow_limits[0]:islow_limits[1], ifast_limits[0]:ifast_limits[1]]
  target_box_size = (islow_limits[1]-islow_limits[0])*(ifast_limits[1]-ifast_limits[0])
  print (-1.0 in all_lunus,"Target=%d"%(A.size))

  minslow=islow_limits[0]
  maxslow=islow_limits[1]
  minfast=ifast_limits[0]
  maxfast=ifast_limits[1]

  B = A!=-1 # elements not flagged by lunus as invalid
  actual_basis_size = np.count_nonzero(B==True)
  if actual_basis_size == target_box_size: return # for now, nothing to do
  cyclic_expansion = -1
  while actual_basis_size < target_box_size: # we will sample from the edges until completing target box
    cyclic_expansion += 1
    cyclic_choice = cyclic_expansion%4
    print("cyclic_choice",cyclic_choice)
    if cyclic_choice == 0:
      if minslow > 0: minslow -= 1
      else: continue
      for ifast in range(minfast,maxfast):
        if self.lunus_filtered_data[ipanel, minslow, ifast]!=-1.0:
          print("Accepted cycle",cyclic_choice,ipanel, minslow, ifast)
          actual_basis_size += 1
          if actual_basis_size==target_box_size: break
    elif cyclic_choice == 1:
      if maxfast + 1 < self.lunus_filtered_data.shape[2]: maxfast += 1 # shape is 2D
      else: continue
      for islow in range(minslow,maxslow):
        if self.lunus_filtered_data[ipanel, islow, maxfast - 1]!=-1.0:
          print("Accepted cycle",cyclic_choice,ipanel, islow, maxfast -1)
          actual_basis_size += 1
          if actual_basis_size==target_box_size: break

#start here continue extending the target box.  Not sure if this code will be sustained

  print ("Success")
  input()
  if (-1.0 in all_lunus): input()
