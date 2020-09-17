from __future__ import absolute_import,print_function, division
import matplotlib.pyplot as plt
import sys,os,time
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

def get_MPI_time():
  with open("rank_0.log","r") as F:
    lines = F.read().strip().split("\n")
    assert "elapsed after Python imports" in lines[-1]
    MPIlog_tokens = lines[-1].split()
    elapsed = float(MPIlog_tokens[11])
    ctime = " ".join(MPIlog_tokens[4:6])
    return time.mktime(time.strptime(ctime,"%Y-%m-%d %H:%M:%S.%f")), elapsed
def get_py_time():
  with open("rank_0.log","r") as F:
    lines = F.read().strip().split("\n")
    assert "elapsed after srun startup" in lines[-2]
    MPIlog_tokens = lines[-2].split()
    elapsed = float(MPIlog_tokens[11])
    ctime = " ".join(MPIlog_tokens[4:6])
    return time.mktime(time.strptime(ctime,"%Y-%m-%d %H:%M:%S.%f")), elapsed
def get_log():
  log_dir = os.path.basename(os.path.abspath("."))
  log_file = os.path.join("..","job%s.out"%(log_dir))
  with open(log_file,"r") as F:
    lines = F.read().strip().split("\n")
    for line in lines:
      if "Started" in line:
        start_time = time.mktime(time.strptime(" ".join(line.split()[2:]),"%a %b %d %H:%M:%S %Y"))
      if "Terminated" in line:
        end_time = time.mktime(time.strptime(" ".join(line.split()[2:]),"%a %b %d %H:%M:%S %Y"))
    print("OK")
  return float(start_time), float(end_time)


def run(params):
  script_start, script_finis = get_log()

  counter = 0
  datum = None
  root=params.input_path
  fig_object = plt.figure()
  good_total = fail_total = 0
  good_timepoints = flex.double()
  good_elapsed = flex.double()
  all_rank = flex.int()

  for filename in os.listdir(root):
    if os.path.splitext(filename)[1] != '.log': continue
    if 'rank' not in filename: continue

    rank = int(filename.split('_')[1].split('.')[0])
    counter += 1
    print (filename, rank)
    for line in open(os.path.join(root,filename)):
      if not line.startswith('idx------finis-------->'): continue
      try:
        _, _, _, _, ts, _, elapsed = line.strip().split()
        epoch_finis = float(ts)
      except ValueError:
        continue
      elapsed = float(elapsed)
      epoch_start = epoch_finis - elapsed
      if datum is None:
        datum = epoch_start
      datum = min(datum, epoch_start)

      good_timepoints.append(epoch_finis)
      good_elapsed.append(elapsed)
      all_rank.append(rank)


  plt.plot(good_timepoints-datum, all_rank, 'g.', markersize="1")
  good_total = len(good_timepoints)
  max_rank = max(all_rank)

  sorted_elapsed = sorted(list(good_elapsed))
  print ("the median weather time is %.5f"%(sorted_elapsed[len(sorted_elapsed)//2]))
  print("Five number summary of %d good image processing times:"%good_total, ["%.5f"%a for a in five_number_summary(good_elapsed)])

  plt.plot([0., max(good_timepoints)-datum],[-(1./30.)*max_rank,-(1./30.)*max_rank], 'r-', label="foreach image")
  print ("The total envelope time is %.1f seconds"%(max(good_timepoints)-datum))
  plt.xlabel('Wall time (sec)')
  plt.ylabel('MPI Rank Number')

  mpi_finish, mpi_elapse = get_MPI_time()
  mpi_start = mpi_finish - mpi_elapse
  plt.plot([mpi_finish-mpi_elapse-datum, mpi_finish-datum],[-(2./30.)*max_rank,-(2./30.)*max_rank], color = "orange", label="MPI comm")
  print ("The total MPI communicator time is %.1f seconds, with %.1f sec before 'foreach' and %.1f sec trailing"%(
           mpi_elapse, datum - mpi_start, mpi_finish - max(good_timepoints) ))


  py_finish, py_elapse = get_py_time()
  py_start = py_finish - py_elapse
  plt.plot([py_finish-py_elapse-datum, py_finish-datum],[-(3./30.)*max_rank,-(3./30.)*max_rank], color = "blue", label="Python time")
  print ("The total Python time is %.1f seconds, with %.1f sec for imports and %.1f sec trailing"%(
           py_elapse, mpi_start - py_start, py_finish - mpi_finish ))

  plt.plot([script_start-datum, script_finis-datum],[-(4./30.)*max_rank,-(4./30.)*max_rank], color = "magenta", label="jsrun time")
  print ("The total script time is %.1f seconds, with %.1f sec for ahead and %.1f sec trailing"%(
           script_finis - script_start, py_start - script_start, script_finis - py_finish ))


  plt.title(params.plot_title)
  if params.pickle_plot:
    from libtbx.easy_pickle import dump
    dump('%s'%params.pickle_filename, fig_object)
  if params.show_plot:
    plt.legend(loc="upper right")
    plt.show()

if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
