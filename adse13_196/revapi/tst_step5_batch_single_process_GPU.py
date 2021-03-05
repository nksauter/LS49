from __future__ import division, print_function
import os
"""Basic test of the exascale API (13_228) applied to the single-panel LS49 step5 simulation.
Test the simulated images against references using the compare script, and
scrape the logs to make a weather plot (not output except in plot mode)."""

if __name__=="__main__":
  os.environ["DEVICES_PER_NODE"] = "1"
  os.environ["N_SIM"] = "10" # total number of images to simulate
  os.environ["USE_EXASCALE_API"] = "True" # "True" or "False" use granular host/device memory transfer
  os.environ["LOG_BY_RANK"] = "1" # Use Aaron's rank logger
  os.environ["RANK_PROFILE"] = "0" # 0 or 1 Use cProfiler, default 1
  os.environ["ADD_SPOTS_ALGORITHM"] = "cuda" # cuda or JH or NKS
  os.environ["ADD_BACKGROUND_ALGORITHM"] = "cuda" # cuda or jh or sort_stable
  os.environ["CACHE_FHKL_ON_GPU"] = "True" # "True" or "False" use single object per rank
  os.environ["MOS_DOM"] = "25"

  # defined script-specific environment variables before application imports
  from LS49.adse13_196.revapi import step5_batch as s5b

  # create and change directory
  os.makedirs("s5bdir",exist_ok=True)
  os.chdir("s5bdir")
  s5b.run_step5_batch(test_without_mpi=True)
  # now compare the results
  from LS49.adse13_196.compare import parse_input,compare_two_images
  params,options = parse_input()
  simulation_template = os.path.join(params.testdir, params.template)
  reference_template = os.path.join(params.refdir, params.template)
  actual_compare = 0
  for i in range(params.ntest):
    newfile = simulation_template%i
    oldfile = reference_template%i

    #print(oldfile,newfile)
    import os
    if not os.path.isfile(oldfile):continue
    if not os.path.isfile(newfile):continue
    print("COMPARISON",i)
    compare_two_images(reference=oldfile, test=newfile, tolerance_delta=params.maxdelta,
      tolerance_count=params.maxcount, verbose_pixels=params.verbose, verbose=False)
    actual_compare += 1
  assert actual_compare >= params.mincompare, "Did not compare %d images, only %d"%(params.mincompare,actual_compare)
  print("OK")
  # finally print OK
