from __future__ import division, print_function
import os
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra

ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data

def run_monochromatic():
  from LS49.sim.step5_pad import tst_all
  tst_all(quick=True)

def compare_two_images(reference, test):
  from dxtbx.format.Registry import Registry
  beam=[]
  data = []
  detector =[]
  headers = []
  for i,fk in enumerate([reference,test]):
    format_instance = Registry.find(fk)
    instance = format_instance(fk)
    beam.append( instance.get_beam() )
    detector.append( instance.get_detector() )
    data.append( instance.get_raw_data() )
    if True: #optional test
      print (beam[-1])
      print (instance.get_goniometer())
      print (detector[-1])
      print (instance.get_scan())
    headers.append( instance.get_smv_header(fk) )

  if  headers[0]==headers[1]:
    print ("Both headers identical")
  else:
    #print headers[0]
    #print headers[1]
    for key in headers[0][1]:
      if key not in headers[1][1]:  print ("second data lacks key",key)
      elif headers[0][1][key] != headers[1][1][key]:
        print ("Key comparison:",key,headers[0][1][key], headers[1][1][key])

  assert len(data[1])==len(data[0])
  diff_data = data[1]-data[0]
  no_differences=True
  for i in diff_data:
    if i!=0:
      print (i,end=" ")
      no_differences=False
  assert no_differences

if __name__=="__main__":
  run_monochromatic()
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000.img.gz"), test="./step5_000000.img.gz")
  print("OK")
