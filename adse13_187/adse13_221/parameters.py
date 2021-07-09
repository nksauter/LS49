from __future__ import division
from dials.array_family import flex
import random, math, numpy as np
from scitbx.matrix import col,sqr
from scipy import constants
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt
import pickle, copy

class variable_mosaicity:
  def __init__(self,value,label,params):
    self.ref_value = value
    self.accepted = flex.int()
    self.running = flex.double()
    self.chain = {label:flex.double()}
    self.proposal = self.ref_value
    CDF_sigma = 1. - math.exp(-0.5) # the CDF at the current position x=sigma
    self.hyperparameter = params.hyperparameter # allowable half width CDF range for the next proposal
    self.target_interval = (CDF_sigma - self.hyperparameter, CDF_sigma + self.hyperparameter)
    self.display_n = 1
    self.display_labels = [label]
    self.formatt = "%s(°)"
    self.label = label
    self.display_ranges = {label:[0, 0.1]} # range in degrees

  def accept(self):
    self.chain[self.label].append(self.proposal)
    self.accepted.append(1)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))
    print("ACCEPTED ",flex.sum(self.accepted),"mosaicity propoasls of ",len(self.accepted))

  def reject(self):
    self.chain[self.label].append(self.chain[self.label][-1])
    self.accepted.append(0)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))

  def generate_next_proposal(self):
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
    self.proposal_vec = [self.proposal] # to provide a uniform interface

  def __del__(self):
    last_half = int(len(self.chain[self.label])//2)
    if last_half<50: return
    stats = flex.mean_and_variance(self.chain[self.label][last_half:])
    print ("""mosaicity %s. initial %f°, final %f±%f°"""%(
    self.label,self.ref_value,stats.mean(), stats.unweighted_sample_standard_deviation(),
     ))

class covariant_cell:
  def __init__(self,ref_crystal, params):
    self.R = ref_crystal
    self.ref_uc = self.R.get_unit_cell().parameters()
    self.accepted = flex.int()
    self.running = flex.double()
    self.transition_probability_ratio = 1.0 # q(X|Y)/q(Y|X), Y=proposal, X=last value
    self.params = params
    self.ki = dict(a=0,b=1,c=2,alpha=3,beta=4,gamma=5)
    self.system = self.R.get_space_group().crystal_system()
    if self.system == "Hexagonal":
      self.proposal_vec = col([self.ref_uc[0],self.ref_uc[2]])
    elif self.system == "Orthorhombic":
      self.proposal_vec = col([self.ref_uc[0],self.ref_uc[1],self.ref_uc[2]])

  @classmethod
  def from_covariance(cls, ref_crystal, params):
    with open(params.covariance,"rb") as M:
      cov = pickle.load(M)
      empcov = cov["populations"].fit_components[0]
      # have the location (a,c) and covariance matrix.
      # want to calculate P(Y|X) and P(X|Y).  Maybe score() gives log-likelihood, if it permits n-sample==1
          # or maybe I need to use the mahalanobis distance in a formula.
      # want to randomly sample a vector from the distribution
    new_instance = cls(ref_crystal, params)
    new_instance.cluster_mean = empcov.location_
    new_instance.cluster_covariance = empcov.covariance_
    new_instance.hyperparameter = params.hyperparameter
    new_instance.display_n = len(cov["features"])
    new_instance.display_labels = cov["features"]
    new_instance.formatt = "%s(Å)"
    if new_instance.system == "Hexagonal":
      assert cov["features"]==['a','c']
    elif new_instance.system == "Orthorhombic":
      assert cov["features"]==['a','b','c']
    new_instance.chain = dict([(L,flex.double()) for L in new_instance.display_labels])
    new_instance.display_ranges = {}
    for idr in range(new_instance.display_n):
      diag_elem = math.sqrt(empcov.covariance_[idr,idr])
      new_instance.display_ranges[new_instance.display_labels[idr]] = [
        new_instance.proposal_vec[idr]-4.*diag_elem, new_instance.proposal_vec[idr]+4*diag_elem]
    return new_instance

  def generate_next_proposal(self):
    MV_vec = col(np.random.multivariate_normal(mean=[0.]*self.display_n, cov=self.cluster_covariance))
    self.proposal_vec = col([self.chain[label][-1] for label in self.display_labels]) + self.hyperparameter * MV_vec
    print(
    "next cell proposal %s"%" ".join(
      ["%s %.6f"%(L, self.proposal_vec[ipv]) for ipv,L in enumerate(self.display_labels)])
    )

  def get_current_crystal_model(self, old_crystal_model):
    from cctbx.uctbx import unit_cell
    if self.system == "Hexagonal":
      old_crystal_model.set_unit_cell(unit_cell((self.proposal_vec[0], self.proposal_vec[0], self.proposal_vec[1], 90., 90., 120.)))
    elif self.system == "Orthorhombic":
      old_crystal_model.set_unit_cell(unit_cell((self.proposal_vec[0], self.proposal_vec[1], self.proposal_vec[2], 90., 90., 90.)))
    else: raise NotImplementedError(self.system)
    return old_crystal_model

  def accept(self):
    for idr in range(self.display_n):
      self.chain[self.display_labels[idr]].append(self.proposal_vec[idr])
    self.accepted.append(1)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))

  def reject(self):
    for idr in range(self.display_n):
      self.chain[self.display_labels[idr]].append(self.chain[self.display_labels[idr]][-1])
    self.accepted.append(0)
    self.running.append(flex.sum(self.accepted)/len(self.accepted))

  def __del__(self):
    last_half = int(len(self.chain[self.display_labels[0]])//2)
    if last_half<50: return
    messages = []
    for idr in range(self.display_n):
      stats = flex.mean_and_variance(self.chain[self.display_labels[idr]][last_half:])
      messages.append("""cell %s.  initial %8.2fÅ, final %8.2f±%4.2fÅ"""%(
        self.display_labels[idr],
        self.ref_uc[self.ki[self.display_labels[idr]]],
        stats.mean(), stats.unweighted_sample_standard_deviation(),
      ))
    print ("\n".join(messages))

class covariant_rot(covariant_cell):
  def __init__(self,ref_crystal, params):
    self.R = copy.deepcopy(ref_crystal) # hold on to a copy of the reference crystal for future
    self.ref_U = sqr(self.R.get_U())
    self.accepted = flex.int()
    self.running = flex.double()
    self.transition_probability_ratio = 1.0 # q(X|Y)/q(Y|X), Y=proposal, X=last value
    self.params = params
    self.ki = dict(x=0,y=1,z=2)

    self.cluster_covariance = np.array([[params.sigmas[0]**2, 0., 0.],
                                        [0., params.sigmas[1]**2, 0.],
                                        [0., 0., params.sigmas[2]**2]])
    self.hyperparameter = params.hyperparameter
    self.display_n = 3
    self.display_labels = ["x","y","z"]
    self.formatt = "rot %s(°)"
    self.chain = dict([(L,flex.double()) for L in self.display_labels])
    self.display_ranges = {}
    for idr in range(self.display_n):
      diag_elem = math.sqrt(self.cluster_covariance[idr,idr])
      self.display_ranges[self.display_labels[idr]] = [
        -2.*diag_elem, +2*diag_elem]
    self.proposal_vec = col(params.value)

  def generate_next_proposal(self):
    MV_vec = col(np.random.multivariate_normal(mean=[0.]*self.display_n, cov=self.cluster_covariance))
    self.proposal_vec = col([self.chain[label][-1] for label in self.display_labels]) + self.hyperparameter * MV_vec
    print(
    "next rotational proposal %s"%" ".join(
      ["%s %8.5f"%(L, self.proposal_vec[ipv]) for ipv,L in enumerate(self.display_labels)])
    )

  def get_current_crystal_model(self, old_crystal_model):
    Rx = col((1.,0.,0.)).axis_and_angle_as_r3_rotation_matrix(angle=self.proposal_vec[0], deg=True)
    Ry = col((0.,1.,0.)).axis_and_angle_as_r3_rotation_matrix(angle=self.proposal_vec[1], deg=True)
    Rz = col((0.,0.,1.)).axis_and_angle_as_r3_rotation_matrix(angle=self.proposal_vec[2], deg=True)
    newU = Rz * (Ry * (Rx * self.ref_U))
    old_crystal_model.set_U(newU)
    return old_crystal_model

  def __del__(self):
    last_half = int(len(self.chain[self.display_labels[0]])//2)
    if last_half<50: return
    messages = []
    for idr in range(self.display_n):
      stats = flex.mean_and_variance(self.chain[self.display_labels[idr]][last_half:])
      messages.append("""rot %s.  initial 0°, final %8.4f±%6.4f°"""%(
        self.display_labels[idr],
        stats.mean(), stats.unweighted_sample_standard_deviation(),
      ))
    print ("\n".join(messages))

class covariant_ncells(covariant_cell):
  def __init__(self, params):
    self.R = copy.deepcopy(params.value) # hold on to a copy of the reference crystal for future
    self.accepted = flex.int()
    self.running = flex.double()
    self.transition_probability_ratio = 1.0 # q(X|Y)/q(Y|X), Y=proposal, X=last value
    self.params = params
    self.ki = dict(a=0,b=1,c=2)
    self.cluster_covariance = np.array([[params.sigmas[0]**2, 0., 0.],
                                        [0., params.sigmas[1]**2, 0.],
                                        [0., 0., params.sigmas[2]**2]])
    self.hyperparameter = params.hyperparameter
    self.display_n = 3
    self.display_labels = ["a","b","c"]
    self.formatt = "ncells %s"
    self.chain = dict([(L,flex.double()) for L in self.display_labels])
    self.display_ranges = {}
    for idr in range(self.display_n):
      diag_elem = math.sqrt(self.cluster_covariance[idr,idr])
      self.display_ranges[self.display_labels[idr]] = [
        0,+100*diag_elem]
    self.proposal_vec = col(params.value)

  def generate_next_proposal(self):
    MV_vec = col(np.random.multivariate_normal(mean=[0.]*self.display_n, cov=self.cluster_covariance))
    print ("THE MC vector",MV_vec)
    print ("The excursion",self.hyperparameter * MV_vec)
    self.proposal_vec = col([self.chain[label][-1] for label in self.display_labels]) + self.hyperparameter * MV_vec
    self.proposal_vec = col(
     (
       abs(self.proposal_vec[0]),abs(self.proposal_vec[1]),abs(self.proposal_vec[2])))
    print(
    "next ncells proposal %s"%" ".join(
      ["%s %d"%(L, self.proposal_vec[ipv]) for ipv,L in enumerate(self.display_labels)])
    )

  def get_current_model(self):
    current = (int(self.proposal_vec[0]), int(self.proposal_vec[1]), int(self.proposal_vec[2]))
    print("CELLS CURRENT",current)
    return current

  def __del__(self):
    last_half = int(len(self.chain[self.display_labels[0]])//2)
    if last_half<50: return
    messages = []
    for idr in range(self.display_n):
      stats = flex.mean_and_variance(self.chain[self.display_labels[idr]][last_half:])
      messages.append("""ncells %s.  initial %3d, final %5.1f±%.1f"""%(
        self.display_labels[idr],self.R[idr],
        stats.mean(), stats.unweighted_sample_standard_deviation(),
      ))
    print ("\n".join(messages))

if __name__=="__main__":
  with open("alt_cryst","rb") as F:
    ref=  pickle.load(F)
  X = covariant_cell.from_covariance(ref)
