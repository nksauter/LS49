from __future__ import absolute_import,print_function, division
import matplotlib.pyplot as plt
import sys,os,time
from libtbx.phil import parse
from libtbx.utils import Sorry

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
  prefix = job
    .type = str
    .help = Use "job" for Summit/LFS, "slurm" for Cori/Slurm
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
  show_finish = False
    .type = bool
    .help = flag to show completion of each task
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
  root = params.input_path
  allrank=[]
  alltime=[]
  allfinrank=[]
  allfintime=[]
  for filename in os.listdir(root):
    if os.path.splitext(filename)[1] != '.log': continue
    if 'rank' not in filename: continue

    rank = int(filename.split('_')[1].split('.')[0])
    for line in open(os.path.join(root,filename)):
        if line.startswith('idx------start'):
          #print (rank,line)
          tokens = line.strip().split()
          stime = float(tokens[4])
          allrank.append(rank)
          alltime.append(stime)
        if line.startswith('idx------finis'):
          tokens = line.strip().split()
          stime = float(tokens[4])
          allfinrank.append(rank)
          allfintime.append(stime)

  begin = min(alltime)
  alltime = [ix -begin for ix in alltime]
  allfintime = [ix -begin for ix in allfintime]

  if params.show_finish: plt.plot(allfintime, allfinrank, "rx")
  plt.plot(alltime, allrank, "g.")
  plt.title("Starting times for %d events in %s"%(len(alltime),root))
  plt.xlabel("Time (s)")
  plt.ylabel("Rank")
  plt.show()

if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:])
  if not params.show_plot:
    import matplotlib
    matplotlib.use("pdf")
  run(params)
