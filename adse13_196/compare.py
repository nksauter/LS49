from __future__ import division
from LS49.tests.tst_monochromatic_image import compare_two_images
import os

def parse_input():
  from iotbx.phil import parse
  master_phil="""
  testdir = .
      .type = path
      .help = the directory with the image simulations to be tested
  refdir = None
      .type = path
      .help = the directory with the reference images
  template = "step5_MPIbatch_%06d.img.gz"
      .type = str
      .help = the file name template
  ntest = 240
      .type = int
      .help = the number of images to test
  maxdelta = 1
      .type = int
      .help = maximum allowable delta value for an image pixel
  maxcount = 100
      .type = int
      .help = maximum number of pixels where a 0 < delta <= maxdelta deviation is acceptable
  verbose = False
      .type = bool
  mincompare = 1
      .type = int
  """
  phil_scope = parse(master_phil)
  # The script usage
  import libtbx.load_env # implicit import
  help_message = '''ADSE13-196.'''
  usage = ""
  '''Initialize the script.'''
  from dials.util.options import ArgumentParser
  # Create the parser
  parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        epilog=help_message)
  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)

  #Program defaults
  if params.refdir is None:
    from LS49 import ls49_big_data
    params.refdir = os.path.join(ls49_big_data,"adse13_196","report1")
  assert os.path.isdir(params.testdir)
  assert os.path.isdir(params.refdir)
  return params,options

if __name__=="__main__":
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
