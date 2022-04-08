from __future__ import absolute_import, division, print_function

import iotbx.phil
from dials.util.options import ArgumentParser, flatten_datablocks
from libtbx.utils import Sorry
from scitbx.array_family import flex

help_message = '''

This program can be used to estimate the gain of the detector. For pixel array
detectors the gain is usually set to 1.00. This means that the pixels behave
according to Poisson statistics. However, for older CCD detectors the gain may
have a different value. This value is important because it can affect, amongst
other things, the ability of the spot finding algorithm which can result in
noise being identified as diffraction spots.

Examples::

  dials.estimate_gain datablock.json

'''

phil_scope = iotbx.phil.parse("""\
  kernel_size = 10,10
    .type = ints(size=2, value_min=1)
  max_images = 1
    .type = int
    .help = "For multi-file images (NeXus for example), report a gain for each"
            "image, up to max_images, and then report an average gain"
  output {
    gain_map = None
      .type = str
      .help = "Name of output gain map file"
  }
""", process_includes=True)

def estimate_gain(imageset, kernel_size=(10,10), output_gain_map=None, max_images = 1):
  detector = imageset.get_detector()

  from dials.algorithms.image.threshold import DispersionThresholdDebug
  gains = flex.double()

  for image_no in xrange(len(imageset)):
    raw_data = imageset.get_raw_data(image_no)
    #from IPython import embed; embed()
    #this_data = raw_data[0]
    #raw_data = (this_data + 80),
    NSQ = 200
    small_section = raw_data[0].matrix_copy_block(400,400,NSQ,NSQ)
    print ("This small section",len(small_section),"mean ist", flex.mean(small_section.as_double()))
    raw_data = (small_section,)


    gain_value = 1
    gain_map = [flex.double(raw_data[i].accessor(), gain_value)
                for i in range(len(detector))]

    mask = imageset.get_mask(image_no)
    mask = (mask[0].matrix_copy_block(400,400,NSQ,NSQ)),
    #from IPython import embed; embed()
    min_local = 0

    # dummy values, shouldn't affect results
    nsigma_b = 6
    nsigma_s = 3
    global_threshold = 0

    kabsch_debug_list = []
    for i_panel in range(len(detector)):
      kabsch_debug_list.append(
        DispersionThresholdDebug(
          raw_data[i_panel].as_double(), mask[i_panel], gain_map[i_panel],
          kernel_size, nsigma_b, nsigma_s, global_threshold, min_local))


    dispersion = flex.double()
    for ipix in range(5,NSQ-15):
      for spix in range(5,NSQ-15):
        data = small_section.matrix_copy_block(ipix,spix,10,10).as_double()
        datasq = data*data
        means = flex.mean(data)
        var = flex.mean(datasq) - (means)**2
        #print(ipix,spix,var,var/means)
        dispersion.append(var/means)

    if True:
     dispersion = flex.double()
     for kabsch in kabsch_debug_list:
      a_section = kabsch.index_of_dispersion().matrix_copy_block(5,5,NSQ-15,NSQ-15)
      print ("mean of a_section",flex.mean(a_section))
      dispersion.extend(a_section.as_1d())

    #ST = flex.mean_and_variance(dispersion)
    #from IPython import embed; embed()

    sorted_dispersion = flex.sorted(dispersion)
    from libtbx.math_utils import nearest_integer as nint

    q1 = sorted_dispersion[nint(len(sorted_dispersion)/4)]
    q2 = sorted_dispersion[nint(len(sorted_dispersion)/2)]
    q3 = sorted_dispersion[nint(len(sorted_dispersion)*3/4)]
    iqr = q3-q1

    print("q1, q2, q3: %.2f, %.2f, %.2f" %(q1, q2, q3))
    if iqr == 0.0:
      raise Sorry('Unable to robustly estimate the variation of pixel values.')

    inlier_sel = (sorted_dispersion > (q1 - 1.5*iqr)) & (sorted_dispersion < (q3 + 1.5*iqr))
    sorted_dispersion = sorted_dispersion.select(inlier_sel)
    gain = sorted_dispersion[nint(len(sorted_dispersion)/2)]
    print("Estimated gain: %.2f" % gain)
    gains.append(gain)

    if image_no == 0:
      gain0 = gain
    if image_no+1 >= max_images:
      break

  if len(gains) > 1:
    stats = flex.mean_and_variance(gains)
    print("Average gain: %.2f +/- %.2f"%(stats.mean(),
      stats.unweighted_sample_standard_deviation()))

  if output_gain_map:
    if len(gains) > 1:
      raw_data = imageset.get_raw_data(0)
    # write the gain map
    import six.moves.cPickle as pickle
    gain_map = flex.double(flex.grid(raw_data[0].all()), gain0)
    with open(output_gain_map, "wb") as fh:
      pickle.dump(gain_map, fh, protocol=pickle.HIGHEST_PROTOCOL)

  if 0:
    sel = flex.random_selection(population_size=len(sorted_dispersion), sample_size=10000)
    sorted_dispersion = sorted_dispersion.select(sel)

    from matplotlib import pyplot
    pyplot.scatter(range(len(sorted_dispersion)), sorted_dispersion)
    pyplot.ylim(0, 10)
    pyplot.show()

  return gain0


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  usage = "%s [options] datablock.json" %libtbx.env.dispatcher_name

  parser = ArgumentParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    check_format=True,
    read_datablocks_from_images=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)

  ## Configure the logging
  #log.config(
    #params.verbosity, info='dials.estimate_gain.log', debug='dials.estimate_gain.debug.log')

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    print('The following parameters have been modified:\n')
    print(diff_phil)

  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    return
  elif len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  assert len(imagesets) == 1
  imageset = imagesets[0]
  estimate_gain(imageset, params.kernel_size, params.output.gain_map, params.max_images)

  return


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
