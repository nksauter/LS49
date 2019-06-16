from __future__ import print_function, division
from six.moves import range
from six.moves import cPickle as pickle
from scitbx.array_family import flex
import scitbx
import math
from LS49.work2_for_aca_lsq.abc_background import fit_roi_multichannel # implicit import
# multichannel needed for unpickling

# %%% boilerplate specialize to packaged big data %%%
import os
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
# %%%%%%
from LS49.sim.step5_pad import data
local_data = data()
Fe_oxidized_model = local_data.get("Fe_oxidized_model")
Fe_reduced_model = local_data.get("Fe_reduced_model")
Fe_metallic_model = local_data.get("Fe_metallic_model")

abc_dir = [ "/net/dials/raid1/sauter/abc_coverage_allrestr",
 #"/net/dials/raid1/sauter/paper1/abc_coverage_superpower_postrefine",
 "/net/dials/raid1/sauter/paper1/abc_coverage_dials_refine",
 "/net/dials/raid1/sauter/paper1/abc_coverage_coarse_ground_truth"]
def get_item(key,trial=0):
      with open("%s/abcX%06d.pickle"%(abc_dir[trial],key),"rb") as F:
        T = pickle.load(F)
      return T

def pprint(M):
  print(M.focus())
  islow,ifast=M.focus()
  for x in range(3,islow-3):
    print (" ".join([("%6.1f"%(M[(x,y)])) for y in range(min(25,ifast))]))

def pprint3(M):
  print(M.focus())
  islow,ifast=M.focus()[1],M.focus()[2]
  for x in range(islow):
    print (" ".join([("%6.1f"%(M[(0,x,y)])) for y in range(min(25,ifast))]))

def show_residual(self,plot=True):
    F = self.sb_data.focus()
    if plot:
      # Residual
      for x in range(3,F[1]-3):
        for y in range(F[2]):
          background = self.a[0]*x+self.a[1]*y+self.a[2]
          model = background + self.a[3] * self.roi[x,y]
          print ("%6.0f"%(self.sb_data[0,x,y]-model),end=' ')
        print()
      print()

from scitbx.array_family.flex import double as flexd
flexd.__str__ = lambda x: str([i for i in x])
from scitbx.array_family.flex import float as flexf
flexf.__str__ = lambda x: str([i for i in x])

def run(myrank):
  T = get_item (myrank)
  for spot in T:
    for item in ['a',
 'asu_idx_C2_setting',
 'bkgrd_a',
 'image_no',
 'n',
 'orig_idx_C2_setting',
 'simtbx_P1_miller',
 'simtbx_intensity_7122',
 'x']:
      print (item, getattr(spot,item))
    for item in ['roi']:
      print (item)
      pprint (getattr(spot, item))
    others = [ "sb_data" ]
    for item in others:
      print (item)
      pprint3 (getattr(spot, item))
    chs = getattr(spot, "channels")
    for c in chs:
      print(c)
      pprint(chs[c])
    #from IPython import embed; embed()

def runs(myrank):
  T = get_item (myrank,0),get_item (myrank,1),get_item (myrank,2)
  for spot in zip(T[0],T[1],T[2]):
    for item in ['a',
 'asu_idx_C2_setting',
 'bkgrd_a',
 'image_no',
 'n',
 'orig_idx_C2_setting',
 'simtbx_P1_miller',
 'simtbx_intensity_7122',
 'x']:
      print (item, getattr(spot[0],item))
      print (item, getattr(spot[1],item))
      print (item, getattr(spot[2],item))

    for item in ['roi']:
      print (item)
      pprint (getattr(spot[0], item))
      print()
      show_residual(spot[0])
      print()
      pprint (getattr(spot[1], item))
      print()
      show_residual(spot[1])
      print()
      pprint (getattr(spot[2], item))
      print()
      show_residual(spot[2])
      print()
    if False:
      others = [ "sb_data" ]
      for item in others:
        print (item)
        pprint3 (getattr(spot[0], item))
        print()
        pprint3 (getattr(spot[1], item))
        print()
        pprint3 (getattr(spot[2], item))
        print()
    continue
    chs = getattr(spot[0], "channels")
    for c in chs:
      print(c)
      pprint(getattr(spot[0], "channels")[c] )
      print()
      pprint(getattr(spot[1], "channels")[c] )
      print()
      pprint( getattr(spot[2], "channels")[c] )
    #from IPython import embed; embed()
"""Results from three trials:
1) Jan 15 "allrestr" dials-refined model
2) Jun 10 dials-refined model
3) Jun 10 ground truth model
Channels:  individual energy-channel models are identical between (1) and (3), suggesting that either the first-draft manuscript
wan't really testing the dials model, or that the June 10 wasn't really ground truth.  June 10 dials model is different intensity
and slightly different position.
ROI:  don't know what it is, but result is similar to channels.
SB_DATA:  Both June 10 models are the same; January different.  Is the raw data?
"""


if __name__=="__main__":
  runs(myrank=12345)
