from __future__ import division, absolute_import, print_function
from libtbx.phil import parse
import libtbx.load_env
from scitbx.matrix import sqr

help_message = '''
libtbx.python test.py input.experiments=idx-step5_000000_integrated_experiments.json
'''

# Create the phil parameters
phil_scope = parse('''
max_hierarchy_level=Auto
  .type = int
  .help = Maximum hierarchy level to compute shifts to
''', process_includes=True)

class Script:
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s experiment1.json experiment2.json reflections1.pickle reflections2.pickle" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_datablocks=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''


    T1="idx-step5_000000_integrated_experiments.json"

    from dxtbx.model.experiment_list import ExperimentListFactory

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    print(params.input.experiments)

    EC = ExperimentListFactory.from_json_file(T1,check_format=False)[0].crystal
    EC.show()
    direct_A = EC.get_A_inverse_as_sqr()
    print(direct_A)
    permute = sqr((0,0,1,0,1,0,-1,0,0))
    sim_compatible = direct_A*permute # permute columns when post multiplying
    print(sim_compatible)

if __name__ == '__main__':

    script = Script()
    script.run()
