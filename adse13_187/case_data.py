from __future__ import division
import os, shutil
from LS49 import ls49_big_data
"""Main idea:  these file names represent the top 75 diffracting events in Run795
of the Oct. 2019 SwissFEL Bernina data collection.  Get the file names with function
get_cori(); transfer them to a git-lfs repo with fetch_to_repo(), and retrieve
from the local copy of repo with retrieve().  Following the links here allows us
to retain the mapping between shot # and top diffracting event # """

def get_cori(N):
  print (N)
  data_file = "/global/cfs/cdirs/m3562/der/run795/top_%d.expt"%N
  print (data_file)
  print (os.path.isfile(data_file), os.path.islink(data_file))
  while os.path.islink(data_file):
    data_file = os.path.join(os.path.dirname(data_file),os.readlink(data_file))
    print(data_file)
  return data_file

def lookup_top(NN=75):
  lookup_cori={}
  lookup_repo={}
  for item in range(NN):
    cori_name = get_cori(item)
    lookup_cori[item] = cori_name
    lookup_repo[item] = os.path.basename(cori_name)
  return lookup_cori, lookup_repo

lookup_repo = \
{0: 'run795_shot648_indexed.expt', 1: 'run795_shot613_indexed.expt', 2: 'run795_shot1193_indexed.expt', 3: 'run795_shot458_indexed.expt', 4: 'run795_shot651_indexed.expt', 5: 'run795_shot1362_indexed.expt', 6: 'run795_shot685_indexed.expt', 7: 'run795_shot1128_indexed.expt', 8: 'run795_shot645_indexed.expt', 9: 'run795_shot447_indexed.expt', 10: 'run795_shot1344_indexed.expt', 11: 'run795_shot745_indexed.expt', 12: 'run795_shot148_indexed.expt', 13: 'run795_shot1349_indexed.expt', 14: 'run795_shot1114_indexed.expt', 15: 'run795_shot654_indexed.expt', 16: 'run795_shot926_indexed.expt', 17: 'run795_shot1107_indexed.expt', 18: 'run795_shot1150_indexed.expt', 19: 'run795_shot595_indexed.expt', 20: 'run795_shot834_indexed.expt', 21: 'run795_shot488_indexed.expt', 22: 'run795_shot178_indexed.expt', 23: 'run795_shot1268_indexed.expt', 24: 'run795_shot308_indexed.expt', 25: 'run795_shot759_indexed.expt', 26: 'run795_shot1286_indexed.expt', 27: 'run795_shot1386_indexed.expt', 28: 'run795_shot760_indexed.expt', 29: 'run795_shot639_indexed.expt', 30: 'run795_shot1054_indexed.expt', 31: 'run795_shot1092_indexed.expt', 32: 'run795_shot1020_indexed.expt', 33: 'run795_shot1258_indexed.expt', 34: 'run795_shot1446_indexed.expt', 35: 'run795_shot1208_indexed.expt', 36: 'run795_shot970_indexed.expt', 37: 'run795_shot710_indexed.expt', 38: 'run795_shot1118_indexed.expt', 39: 'run795_shot1158_indexed.expt', 40: 'run795_shot1201_indexed.expt', 41: 'run795_shot702_indexed.expt', 42: 'run795_shot598_indexed.expt', 43: 'run795_shot764_indexed.expt', 44: 'run795_shot1301_indexed.expt', 45: 'run795_shot396_indexed.expt', 46: 'run795_shot479_indexed.expt', 47: 'run795_shot684_indexed.expt', 48: 'run795_shot751_indexed.expt', 49: 'run795_shot722_indexed.expt', 50: 'run795_shot1212_indexed.expt', 51: 'run795_shot1322_indexed.expt', 52: 'run795_shot295_indexed.expt', 53: 'run795_shot1214_indexed.expt', 54: 'run795_shot1161_indexed.expt', 55: 'run795_shot12_indexed.expt', 56: 'run795_shot826_indexed.expt', 57: 'run795_shot441_indexed.expt', 58: 'run795_shot601_indexed.expt', 59: 'run795_shot840_indexed.expt', 60: 'run795_shot704_indexed.expt', 61: 'run795_shot699_indexed.expt', 62: 'run795_shot375_indexed.expt', 63: 'run795_shot428_indexed.expt', 64: 'run795_shot69_indexed.expt', 65: 'run795_shot1275_indexed.expt', 66: 'run795_shot712_indexed.expt', 67: 'run795_shot132_indexed.expt', 68: 'run795_shot652_indexed.expt', 69: 'run795_shot1288_indexed.expt', 70: 'run795_shot491_indexed.expt', 71: 'run795_shot837_indexed.expt', 72: 'run795_shot1289_indexed.expt', 73: 'run795_shot28_indexed.expt', 74: 'run795_shot392_indexed.expt'}

def fetch_to_repo(lookup_cori, lookup_repo):
  targetdir = os.path.join(ls49_big_data,"adse13_228")
  for item in lookup_cori:
    print(item)
    _src = lookup_cori[item]
    _dst = os.path.join(targetdir,lookup_repo[item])
    shutil.copy2(src=_src, dst=_dst)

def retrieve_from_repo(N):
  targetdir = os.path.join(ls49_big_data,"adse13_228")
  return os.path.join(targetdir, lookup_repo[N])

if __name__=="__main__":
  lookup_cori, lookup_repo = lookup_top(75)
  #print (lookup_repo)
  fetch_to_repo(lookup_cori, lookup_repo)
