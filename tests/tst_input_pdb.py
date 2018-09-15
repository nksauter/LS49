from __future__ import division, print_function
import os

ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment

def get_results():
  pdb_lines = open(os.path.join(ls49_big_data,"1m2a.pdb"),"r").read()
  verbose=False
  if verbose:
    for line in pdb_lines.split("\n"):
      print(line)
  return pdb_lines

def create_reference_results_from_pdb():
  from mmtbx.command_line.fetch_pdb import run2
  run2(["1m2a"])
  pdb_lines = open("1m2a.pdb").read()
  return pdb_lines


if __name__=="__main__":
  data_store = get_results()
  pdb_refern = create_reference_results_from_pdb()
  assert data_store==pdb_refern
  print("pdb file #characters=",len(data_store))
  print("OK")
