from __future__ import division, print_function

from iotbx.phil import parse

help_message = '''
implement Max Likelihood refinement of the fp, fdp parameters.
'''

master_phil="""
N_total = 2000
    .type = int(value_min=0, value_max=100000)
    .help = How many elements of the total simulation to use for parameter refinement
    .help = Possible max working size for 64-core rack server is 6400
tester{
  rank = 45
    .type = int
    .help = For the purpose of testing a single worker rank, which rank number to simulate
  size = 1024
    .type = int
    .help = For the purpose of testing a single worker rank, what total MPI size to simulate
}
starting_model{
  algorithm = to_file from_file_static from_file_all
  .type = choice
  .help = Compute the energy-dependent intensities and derivatives using chosen procedure
  .help = from_file_static choice indicates the intensities will be computed on the fly for
  .help =   each iteration, however the static Fcalcs still come from the file.
  .help = to_file choice computes the array in rank 0, then pickles it and exits the program
  .help = from_file reads intensities from file for the initial iteration, then compute on the fly
  filename = new_global_fdp_big_data.pickle
      .type = path
      .help = write out
}
"""

phil_scope = parse(master_phil)
