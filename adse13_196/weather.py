from __future__ import absolute_import,print_function, division
import matplotlib.pyplot as plt
import sys,os
from libtbx.phil import parse
from libtbx.utils import Sorry
from scitbx.array_family import flex
from scitbx.math import five_number_summary

message = ''' script to get a sense of the computational performance of every rank while processing data.
              End product is a plot of wall time vs MPI rank number with every data point being that of a frame
              processed by step5_batch.py. The information is read in from the debug files created by
              step5_batch.py.
              Example usage:
              libtbx.python weather.py
'''
phil_scope = parse('''
  input_path = .
    .type = str
    .help = path to where the processing results are. For example path to XXX_rgYYYY
  num_nodes = 1
    .type = int
    .help = Number of nodes used to do data processing. Used in timing information
  num_cores_per_node = 72
    .type = int
    .help = Number of cores per node in the machine (default is for Cori KNL)
  wall_time = 3600
    .type = int
    .help = total wall time (seconds) taken for job to finish. Used for plotting node-partitioning
  plot_title = Computational weather plot
    .type = str
    .help = title of the computational weather plot
  show_plot = True
    .type = bool
    .help = flag to indicate if plot should be displayed on screen
  pickle_plot = False
    .type = bool
    .help = If True, will pickle matplotlib session so that it can be opened later for analysis/viewing \
            https://stackoverflow.com/questions/29160177/matplotlib-save-file-to-be-reedited-later-in-ipython
  pickle_filename = fig_object.pickle
    .type = str
    .help = Default name of pickled matplotlib plot saved to disk
''')

def params_from_phil(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse(file_name=arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  return params

def run(params):
  counter = 0
  reference = None
  root=params.input_path
  fig_object = plt.figure()
  good_total = fail_total = 0
  for filename in os.listdir(root):
    if os.path.splitext(filename)[1] != '.log': continue
    if 'rank' not in filename: continue
    fail_timepoints = []
    good_timepoints = []
    rank = int(filename.split('_')[1].split('.')[0])
    counter += 1
    print (filename, rank)
    for line in open(os.path.join(root,filename)):
      if not line.startswith('idx------finis-------->'): continue
      try:
        _, _, _, _, ts, _, elapsed = line.strip().split()
        ts = float(ts)
      except ValueError:
        continue
      if reference is None:
        reference = ts - float(elapsed)

      status = 'done'
      if status in ['stop','done','fail']:
        if status == 'done':
          good_timepoints.append(ts-reference)
        else:
          fail_timepoints.append(ts-reference)
        ok = True
      else:
        ok = False
    plt.plot(fail_timepoints, [rank]*len(fail_timepoints), 'b.')
    plt.plot(good_timepoints, [rank]*len(good_timepoints), 'g.')
    fail_total += len(fail_timepoints)
    good_total += len(good_timepoints)
    if not ok:
      plt.plot([ts - reference], [rank], 'rx')
    #if counter > 100: break

  fail_deltas = [fail_timepoints[i+1] - fail_timepoints[i] for i in range(len(fail_timepoints)-1)]
  good_deltas = [good_timepoints[i+1] - good_timepoints[i] for i in range(len(good_timepoints)-1)]
  if fail_deltas: print("Five number summary of %d fail image processing times:"%fail_total, five_number_summary(flex.double(fail_deltas)))
  if good_deltas: print("Five number summary of %d good image processing times:"%good_total, five_number_summary(flex.double(good_deltas)))

  for i in range(params.num_nodes):
    plt.plot([0,params.wall_time], [i*params.num_cores_per_node-0.5, i*params.num_cores_per_node-0.5], 'r-')
  plt.xlabel('Wall time (sec)')
  plt.ylabel('MPI Rank Number')
  plt.title(params.plot_title)
  if params.pickle_plot:
    from libtbx.easy_pickle import dump
    dump('%s'%params.pickle_filename, fig_object)
  if params.show_plot:
    plt.show()

if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
