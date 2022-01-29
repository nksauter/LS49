from __future__ import division
import numpy as np
from scitbx.array_family import flex
import scitbx.lbfgs

def is_valid(value): return value>-1.0
def replacement_pixels(lunus_filtered_data, ipanel, islow_limits, ifast_limits, shoebox):
  """main idea: LUNUS marks invalid model pixels as -1.0.  If so, try to find an equal number of valid
  pixels from around the edges of the box, prior to fitting the 2nd-order Taylor function.
  """
  A = all_lunus = lunus_filtered_data[ipanel, islow_limits[0]:islow_limits[1], ifast_limits[0]:ifast_limits[1]]
  target_box_size = (islow_limits[1]-islow_limits[0])*(ifast_limits[1]-ifast_limits[0])
  #print (-1.0 in all_lunus,"Target=%d"%(A.size))
  minslow=islow_limits[0]
  maxslow=islow_limits[1]
  minfast=ifast_limits[0]
  maxfast=ifast_limits[1]

  # Get the initialization data for fitting the pixels on one shoebox
  S = shoebox
  peakslow = int(S.centroid_all().px.position[1]) # approx slow-position of peak
  peakfast = int(S.centroid_all().px.position[0]) # approx fast-position of peak
  FIT = fit_background_one_spot(boxsize = target_box_size, peakslow=peakslow, peakfast=peakfast, gain=9.5)
  for islow in range(islow_limits[0], islow_limits[1]):
    for ifast in range(ifast_limits[0], ifast_limits[1]):
      value = lunus_filtered_data[ipanel,islow,ifast]
      if is_valid(value): FIT.addpixel(islow,ifast,value)

  actual_basis_size = len(FIT.SP) # np.count_nonzero(B==True)
  cyclic_expansion = -1
  while actual_basis_size < target_box_size: # we will sample from the edges until completing target box
    cyclic_expansion += 1
    cyclic_choice = cyclic_expansion%4
    if cyclic_choice == 0:
      if minslow > 0: minslow -= 1
      else: continue
      for ifast in range(minfast,maxfast):
        if is_valid(lunus_filtered_data[ipanel, minslow, ifast]):
          FIT.addpixel(minslow,ifast,value=lunus_filtered_data[ipanel, minslow, ifast])
          #print("Accepted cycle",cyclic_choice,ipanel, minslow, ifast)
          actual_basis_size += 1
          if actual_basis_size==target_box_size: break
    elif cyclic_choice == 1:
      if maxfast + 1 < lunus_filtered_data.shape[2]: maxfast += 1
      else: continue
      for islow in range(minslow,maxslow):
        if is_valid(lunus_filtered_data[ipanel, islow, maxfast - 1]):
          FIT.addpixel(islow,maxfast-1,value=lunus_filtered_data[ipanel, islow, maxfast-1])
          #print("Accepted cycle",cyclic_choice,ipanel, islow, maxfast -1)
          actual_basis_size += 1
          if actual_basis_size==target_box_size: break
    elif cyclic_choice == 2:
      if maxslow + 1 < lunus_filtered_data.shape[1]: maxslow += 1
      else: continue
      for ifast in range(maxfast-1,minfast-1,-1):
        if is_valid(lunus_filtered_data[ipanel, maxslow - 1, ifast]):
          FIT.addpixel(maxslow-1,ifast,value=lunus_filtered_data[ipanel, maxslow-1, ifast])
          #print("Accepted cycle",cyclic_choice,ipanel, maxslow - 1, ifast)
          actual_basis_size += 1
          if actual_basis_size==target_box_size: break
    elif cyclic_choice == 3:
      if minfast > 0: minfast -= 1
      else: continue
      for islow in range(maxslow-1,minslow-1,-1):
        if is_valid(lunus_filtered_data[ipanel, islow, minfast]):
          FIT.addpixel(islow,minfast,value=lunus_filtered_data[ipanel, islow, minfast])
          #print("Accepted cycle",cyclic_choice,ipanel, islow, minfast)
          actual_basis_size += 1
          if actual_basis_size==target_box_size: break

  FIT.initialize_A()
  FIT.initialize_and_run()
  return FIT

class fit_background_one_spot:
  def __init__(self, boxsize, peakslow, peakfast, gain):
    self.n = 6
    self.peakslow = peakslow
    self.peakfast = peakfast
    self.gain = gain
    self.SP = flex.double()
    self.FP = flex.double()
    self.SS = flex.double()
    self.SF = flex.double()
    self.FF = flex.double()
    self.KI = flex.double()
    self.LI = flex.double()
    self.B = 0.; self.C = 0.; self.D = 0.; self.E = 0.; self.F = 0.;
  def addpixel(self, slow, fast, value):
    self.SP.append(slow-self.peakslow)
    self.FP.append(fast-self.peakfast)
    self.SS.append(0.5 * self.SP[-1]**2)
    self.SF.append(0.5 * self.SP[-1] * self.FP[-1])
    self.FF.append(0.5 * self.FP[-1]**2)
    self.KI.append(value/self.gain)
  def initialize_A(self):
    # initial fit is to set the average
    self.A = self.gain * flex.mean(self.KI)
  def model_T(self,slow,fast):
    sprime = slow-self.peakslow; fprime = fast-self.peakfast
    return self.A + sprime*self.B + fprime*self.C + 0.5*(
           sprime*sprime*self.D + sprime*fprime*self.E + fprime*fprime*self.F)
  def initialize_and_run(self):
    self.x = flex.double([self.A,self.B,self.C,self.D,self.E,self.F,])
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-3, #significantly (4x) quicker than 1.e-4
        max_calls=100)
    )
    self.a = self.x
    self.A,self.B,self.C,self.D,self.E,self.F = tuple(self.a)

  def compute_functional_and_gradients(self):
    self.a = self.x
    f = 0.;
    g = flex.double(self.n)
    vector_T = flex.double(len(self.SP),self.x[0]) + self.SP*self.x[1] + self.FP*self.x[2] + 0.5*(
           self.SS*self.x[3] + self.SF*self.x[4] + self.FF*self.x[5])
    vector_lambda = vector_T/self.gain
    if (vector_lambda <= 0).count(True)>0:
      raise RuntimeError("raising exception to avoid log(value<=0)")
    f = flex.sum(vector_lambda - (self.KI * flex.log(vector_lambda)))
    inner_paren = flex.double(len(self.SP),1.) - (self.KI/vector_lambda)
    g_list = [flex.sum( deriv * inner_paren ) for deriv in
               [flex.double(len(self.SP),1.), self.SP, self.FP, self.SS, self.SF, self.FF]]
    #self.print_step("LBFGS stp",f)
    g_list[3]=0.; g_list[4]=0.; g_list[5]=0. # turn off the 2nd-order Taylor term
    g = flex.double(g_list)/self.gain
    return f,g
  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%10.4f"%a for a in self.x]),"]")
