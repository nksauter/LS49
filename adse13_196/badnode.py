from __future__ import absolute_import,print_function, division
import sys,os,time
from libtbx.phil import parse
from libtbx.utils import Sorry
from scitbx.array_family import flex
from scitbx.math import five_number_summary

message = ''' script to get a sense of the computational performance of every rank while processing data.
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
  devices_per_node = 6
    .type = int
  ranks_per_device = 7
    .type = int
  cores_per_node = 42
    .type = int
  ranks = None
    .type = int
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

def get_log():
  log_dir = os.path.basename(os.path.abspath("."))
  log_file = os.path.join("..","job%s.out"%(log_dir))
  with open(log_file,"r") as F:
    lines = F.read().strip().split("\n")
    for line in lines:
      if "hello from rank 1" in line:
        tokens = line.split()[12::2]
        break
    print("OK")
  return tokens

def run(params):
  node_names = get_log()
  counter = 0
  root=params.input_path
  good_total = fail_total = 0
  good_elapsed = flex.double()
  good_channels = flex.double()
  good_logger = flex.double()
  all_rank = flex.int()
  channels_rank = flex.int()
  logger_rank = flex.int()
  device_elapsed = dict()
  for rank in range(params.ranks):
    filename = "rank_%d.log"%(rank)
    if 'rank' not in filename: continue
    node_num = rank//params.ranks_per_device//params.devices_per_node
    device_num = (rank//params.ranks_per_device)%params.devices_per_node
    device_addr = (node_num,device_num)
    device_elapsed[device_addr] = device_elapsed.get(device_addr, [])

    counter += 1
    if counter%100==1: print (filename, counter)
    for line in open(os.path.join(root,filename)):
      if not line.startswith('idx------finis-------->'): continue
      try:
        _, _, _, _, ts, _, elapsed = line.strip().split()
        epoch_finis = float(ts)
      except ValueError:
        continue
      elapsed = float(elapsed)
      device_elapsed[device_addr].append(elapsed)

      good_elapsed.append(elapsed)
      all_rank.append(rank)
    print("Rank",rank,"node",node_num,node_names[node_num],"device",device_num)
  print("There are %d images"%(len(all_rank)))
  print("There are %d device addresses"%(len(device_elapsed)))
  for node_num,device_num in device_elapsed:
    good_elapsed = device_elapsed[(node_num,device_num)]
    sorted_elapsed = sorted(good_elapsed)
    print ("Median elapsed","node",node_num,node_names[node_num],"device",device_num,"is %.4f"%(sorted_elapsed[len(sorted_elapsed)//2]),
           "5# summary of %d times:"%(len(good_elapsed)), ["%.4f"%a for a in five_number_summary(good_elapsed)])

if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
