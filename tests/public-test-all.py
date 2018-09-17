# LIBTBX_SET_DISPATCHER_NAME cctbx_regression.test_all

from __future__ import division
from libtbx.command_line import run_tests_parallel
import sys
import libtbx.load_env # implicit import

if (__name__ == "__main__") :
  args = [
    "module=libtbx",
    "module=boost_adaptbx",
    "module=scitbx",
    "module=cctbx",
    "module=iotbx",
    "module=smtbx",
    "module=spotfinder",
    "module=rstbx",
    "module=fable",
    "module=simtbx",
    "module=dxtbx",
    "module=dials",
    "module=annlib_adaptbx",
    "module=cbflib_adaptbx",
    "module=LS49",
    "nproc=Auto",
  ]
# Also add these if available
#    "module=labelit",
#    "module=labelit_regression",
#    "module=dials_regression",
#    "module=xfel_regression",

  if (run_tests_parallel.run(args) > 0) :
    sys.exit(1)
