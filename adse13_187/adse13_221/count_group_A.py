from __future__ import division
import os
from dials.array_family import flex
dirpath = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b"
from dxtbx.model.experiment_list import ExperimentListFactory

def from_files():
  for x in range(7534):
    if x%500==0: print(x)
    expt_file = "split_%04d.expt"%x
    expts = ExperimentListFactory.from_json_file(os.path.join(dirpath,expt_file),
                                              check_format=False)
    refl_file = "split_%04d.refl"%x
    refl_table = flex.reflection_table.from_file(os.path.join(dirpath,refl_file))

    yield refl_table,expts[0]

def run_job():
  from cctbx.crystal import symmetry
  from cctbx.uctbx import unit_cell
  from cctbx.sgtbx import space_group
  sym = symmetry(unit_cell=unit_cell((78.68, 78.68, 265.33, 90, 90, 120)),
                 space_group="P6522")
  millers = flex.miller_index()
  for ref,exp in from_files():
    idx = ref["miller_index"]
    millers.extend(idx)
  from cctbx.miller import set as mset, array
  MM = mset(crystal_symmetry=sym,anomalous_flag=True,indices=millers) # make miller set from all obs
  MA = MM.map_to_asu() # map to asu
  #dval = MA.unit_cell().d(MA.indices())
  #order = flex.sort_permutation(dval, reverse=True)
  #MAsort = MA.select(order)
  #for item in MAsort.indices(): print(item)
  #print (MAsort.indices().size(), len(set(MAsort.indices())))
  #ZZ = MAsort.array(data=flex.double(MAsort.size(),1.))
  MAobs = MA.array(data=flex.double(MA.size(),1.))
  ME = MAobs.merge_equivalents()
  new_red = ME.redundancies() # all observations with multiplicity of measurement
  z = new_red.complete_array(d_max=100000.,new_data_value=0) # fill in missing observations
  dval = z.unit_cell().d(z.indices())
  order = flex.sort_permutation(dval, reverse=True)
  MAsort = z.select(order)
  print("d_min is %.3f"%(flex.min(MA.d_spacings().data())))
  print("there are altogether %d measurements"%MA.indices().size())
  print("this is uniqueness : %d measurements"%new_red.indices().size())
  print("total in asymmetric unit %d"%MAsort.indices().size())
  print("maximum multiplicity is %d"%(flex.max(MAsort.data())))
  print("the mean multiplicity is %.2f"%(flex.mean(MAsort.data().as_double())))
  zeroes = (MAsort.data()==0).count(True)
  print("unmeasured %d or %.2f%% of total"%(zeroes, 100.*zeroes/MAsort.indices().size()))

  MAd = MA.select(MA.d_spacings().data()>1.63)
  new_redd = new_red.select(new_red.d_spacings().data()>1.63)
  MAsortd = MAsort.select(MAsort.d_spacings().data()>1.63)
  print("Cut, there are altogether %d measurements"%MAd.indices().size())
  print("Cut, this is uniqueness : %d measurements"%new_redd.indices().size())
  print("Cut, total in asymmetric unit %d"%MAsortd.indices().size())
  print("Cut, maximum multiplicity is %d"%(flex.max(MAsortd.data())))
  print("Cut, the mean multiplicity is %.2f"%(flex.mean(MAsortd.data().as_double())))
  zeroes = (MAsortd.data()==0).count(True)
  print("Cut, unmeasured %d or %.2f%% of total"%(zeroes, 100.*zeroes/MAsortd.indices().size()))

  if True:
      from matplotlib import pyplot as plt
      plt.plot(range(MAsort.indices().size()),MAsort.data(),"r.")
      plt.plot(range(MAsortd.indices().size()),MAsortd.data().as_double()-0.1,"b.")
      plt.title("Redundant measurement vs. Bragg spot position")
      plt.xlabel("Spots ordered by increasing Bragg angle →")
      plt.ylabel("Multiplicity of strong spot obs (Group A)")
      plt.show()


if __name__=="__main__":
  run_job()

