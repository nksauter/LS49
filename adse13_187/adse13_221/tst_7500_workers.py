from __future__ import division
from LS49.adse13_187.adse13_221.lunus_wrap import run, parse_input

def top_75_iterator():
  import re
  from LS49.adse13_187.case_data import lookup_repo
  for ikey, key in enumerate(lookup_repo):
    assert ikey==key
    item = lookup_repo[key]
    match = re.search("shot([0-9]*)",item)
    event_idx = int(match.groups()[0])
    yield event_idx

def multiple_cases():
  def get_any_case():
    params,options = parse_input()
    from LS49.adse13_187.adse13_221.ad_hoc_run795_lookup import conversion
    for idx,item in enumerate(top_75_iterator()):
      print("internal",idx,item)
      #if idx>3: exit() # quick check the first one
      run_no = 795 # only look at run number 795
      params.trusted_mask = "/global/cscratch1/sd/nksauter/adse13_187/bernina/trusted_Py3.mask"
      params.refl = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_%04d.refl"%(
        conversion[run_no][item]
        )
      params.expt = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split_cr/split_%04d.expt"%item
      params.output.label = "L3"
      params.output.index = int(idx)
      yield params
  for prm in get_any_case():
    print(prm.expt)
    print(prm.refl)
    print()

def single_case():
  def get_case_1():
      params,options = parse_input()
      params.trusted_mask = "/global/cscratch1/sd/nksauter/adse13_187/bernina/trusted_Py3.mask"
      params.refl = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_0309.refl"
      params.expt = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split_c/split_0648.expt"
      params.output.label = "L4"
      params.output.index = 648
      return params
  prm = get_case_1()
  run(prm)

if __name__ == "__main__":
  # single_case()
  multiple_cases()
