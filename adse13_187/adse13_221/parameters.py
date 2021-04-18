from __future__ import division
from dials.array_family import flex
import random, math, numpy as np
from scitbx.matrix import col
from scipy import constants
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt
import pickle

class variable_mosaicity:
  def __init__(self,value):
    self.ref_value = value
    self.accepted = flex.int()
    self.running = flex.double()
    self.chain= flex.double()
    self.proposal = self.ref_value
    CDF_sigma = 1. - math.exp(-0.5) # the CDF at the current position x=sigma
    self.hyperparameter = 0.2 # allowable half width CDF range for the next proposal
    self.target_interval = (CDF_sigma - self.hyperparameter, CDF_sigma + self.hyperparameter)

  def accept(self):
    self.chain.append(self.proposal)
    self.accepted.append(1)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))
    print("ACCTEPTED ",flex.sum(self.accepted),"mosaicity propoasls of ",len(self.accepted))

  def reject(self):
    self.chain.append(self.chain[-1])
    self.accepted.append(0)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))

  def generate_next_proposal(self):
    #import numpy.random
    #self.proposal = numpy.random.rayleigh(scale=self.chain[-1])
    """let's go back to basics:"""
    last = sigma = float(self.proposal) # current value
    deviate = random.random()
    selected_targetCDF = self.target_interval[0] + deviate * (self.target_interval[1]-self.target_interval[0])
    proposal = self.proposal = math.sqrt(-2.* sigma *sigma * math.log(1.-selected_targetCDF))
    # gives a random Rayleigh deviate in the target interval around present value "sigma"

    # Rayleigh PDF is (x/ss)*exp(-xx/(2ss))
    prob_proposal_given_last = (proposal/(last*last))*math.exp(-1.*proposal*proposal/(2*last*last))
    prob_last_given_proposal = (last/(proposal*proposal))*math.exp(-1.*last*last/(2*proposal*proposal))
    self.transition_probability_ratio = prob_last_given_proposal/prob_proposal_given_last
    # q(X|Y)/q(Y|X), Y=proposal, X=last value

    print('next mosaicity proposal %.6f'%(self.proposal))

  def __del__(self):
    last_half = int(len(self.chain)//2)
    stats = flex.mean_and_variance(self.chain[last_half:])
    print ("""mosaicity 𝜂. initial %f°, final %f±%f°
"""%(self.ref_value,stats.mean(), stats.unweighted_sample_standard_deviation(),
     ))

class variable_cell:
  def __init__(self,ref_crystal):
    self.R = ref_crystal
    """from covariance 78.68±0.04  265.33±0.37
       from alt cell: a=78.598, c= 265.317"""
    self.ref_uc = self.R.get_unit_cell().parameters()
    self.accepted = flex.int()
    self.running = flex.double()
    self.a_chain= flex.double()
    self.a_proposal = self.ref_uc[0]

    self.c_chain= flex.double()
    self.c_proposal = self.ref_uc[2]

    self.a_sigma = 0.04 # 0.01 # 0.04
    self.c_sigma = 0.37 # 0.03 # 0.37
    self.transition_probability_ratio = 1.0 # q(X|Y)/q(Y|X), Y=proposal, X=last value

  def get_current_crystal_model(self):
    from cctbx.uctbx import unit_cell
    self.R.set_unit_cell(unit_cell((self.a_proposal, self.a_proposal, self.c_proposal, 90., 90., 120.)))
    return self.R

  def accept(self):
    self.a_chain.append(self.a_proposal); self.c_chain.append(self.c_proposal)
    self.accepted.append(1)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))

  def reject(self):
    self.a_chain.append(self.a_chain[-1]); self.c_chain.append(self.c_chain[-1])
    self.accepted.append(0)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))

  def generate_next_proposal(self):
    self.a_proposal = random.gauss(mu=self.a_chain[-1], sigma=self.a_sigma)
    self.c_proposal = random.gauss(mu=self.c_chain[-1], sigma=self.c_sigma)
    print('next cell a-c proposal %.6f %.6f'%(self.a_proposal, self.c_proposal))

  def __del__(self):
    last_half = int(len(self.a_chain)//2)
    stats_a = flex.mean_and_variance(self.a_chain[last_half:])
    stats_c = flex.mean_and_variance(self.c_chain[last_half:])
    print ("""cell a.  initial %8.2fÅ, final %8.2f±%4.2fÅ
cell c.  initial %8.2fÅ, final %8.2f±%4.2fÅ
"""%(self.ref_uc[0],stats_a.mean(), stats_a.unweighted_sample_standard_deviation(),
     self.ref_uc[2],stats_c.mean(), stats_c.unweighted_sample_standard_deviation(),
     ))

class covariant_cell (variable_cell):

  @classmethod
  def from_covariance(cls, ref_crystal):
    with open("../covariance_cytochrome.pickle","rb") as M:
      cov = pickle.load(M)
      empcov = cov["populations"].fit_components[0]
      # have the location (a,c) and covariance matrix.
      # want to calculate P(Y|X) and P(X|Y).  Maybe score() gives log-likelihood, if it permits n-sample==1
          # or maybe I need to use the mahalanobis distance in a formula.
      # want to randomly sample a vector from the distribution
    new_instance = cls(ref_crystal)
    new_instance.cluster_mean = empcov.location_
    new_instance.cluster_covariance = empcov.covariance_
    new_instance.hyperparameter = 0.5
    return new_instance

  def generate_next_proposal(self):

    MV_vec = col(np.random.multivariate_normal(mean=[0.,0.], cov=self.cluster_covariance))
    proposal_vec = col([self.a_chain[-1],self.c_chain[-1]]) + self.hyperparameter * MV_vec
    self.a_proposal = proposal_vec[0]
    self.c_proposal = proposal_vec[1]
    print('next cell a-c proposal %.6f %.6f'%(self.a_proposal, self.c_proposal))

if __name__=="__main__":
  with open("alt_cryst","rb") as F:
    ref=  pickle.load(F)
  X = covariant_cell.from_covariance(ref)
