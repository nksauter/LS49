from __future__ import division
from dials.array_family import flex
import os, pickle
from dxtbx.model.experiment_list import ExperimentListFactory
from LS49.adse13_187.adse13_221.data_slots import application_slots

class mask_manager:
  def __init__(self,trusted_mask,refl_table,expt):
    self.trusted_mask = trusted_mask
    self.refl_table = refl_table
    self.expt = expt
    self.view = application_slots
    self.n_panels = len(self.expt.detector)

  @classmethod
  def from_files(cls, trusted_mask_file, refl_file, expt_file):
    with open(trusted_mask_file,"rb") as M:
      mask = pickle.load(M)
    refl_table = flex.reflection_table.from_file(refl_file)
    expts = ExperimentListFactory.from_json_file(expt_file,
                                              check_format=True)
    return cls(mask,refl_table,expts[0])

  @classmethod
  def from_mask_file_expt_refl(cls, trusted_mask_file, expt, refl):
    with open(trusted_mask_file,"rb") as M:
      mask = pickle.load(M)
    refl_table = refl
    return cls(mask,refl_table,expt)

  def get_trusted_and_refl_mask(self):
    mask = self.trusted_mask
    Z = self.refl_table
    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    size = len(Z)
    BORDER=3
    self.shoebox_mask=[]
    for midx in range(self.n_panels):
      self.shoebox_mask.append(flex.bool(flex.grid(mask[midx].focus()), False))
      # should be a False field

    N_true = 0
    N_bad_pizel = 0
    for sidx in range(size): #loop through the shoeboxes
      ipanel = P[sidx]
      false_field = self.shoebox_mask[ipanel]
      slow_size = false_field.focus()[0]
      fast_size = false_field.focus()[1]
      bbox = S[sidx].bbox
      for islow in range(max(0,bbox[2]-BORDER), min(slow_size,bbox[3]+BORDER)):
        for ifast in range(max(0,bbox[0]-BORDER), min(fast_size,bbox[1]+BORDER)):
          false_field[islow*slow_size + ifast]=True
          N_true += 1
          if mask[ipanel][islow*slow_size + ifast] is False: N_bad_pizel+=1
    self.resultant = []
    for midx in range(self.n_panels):
      newmask = mask[midx].__and__(self.shoebox_mask[midx])
      self.resultant.append(newmask)
    print (N_true,"pixels were visited in the %d shoeboxes (with borders)"%size)
    print (N_bad_pizel,"of these were bad pixels")
    self.monolithic_mask_whole_detector_as_1D_bool = flex.bool()
    for ipnl in range(len(self.resultant)):
      pnl = self.resultant[ipnl]
      self.monolithic_mask_whole_detector_as_1D_bool.extend(pnl.as_1d())
    assert len(self.monolithic_mask_whole_detector_as_1D_bool)==self.n_panels*254*254

  def resultant_mask_to_file(self,file_name):
    all_done = tuple(self.resultant)
    with open(file_name,"wb") as M:
      pickle.dump(all_done,M)

  def get_lunus_repl(self):
    # first get the lunus image
    from lunus.command_line.filter_peaks import get_image_params
    from lunus import LunusDIFFIMAGE

    imageset = self.expt.imageset
    data = imageset[0]
    assert isinstance(data, tuple) # assume a tuple of flex::double over detector panels
    # Instantiate a LUNUS diffraction image
    A = LunusDIFFIMAGE(len(data))

    # Populate the image with multipanel data

    for pidx in range(len(data)):
      A.set_image(pidx,data[pidx])

    # Define the LUNUS image parameters
    deck = '''
#lunus input deck
#punchim_xmin=1203
#punchim_ymin=1250
#punchim_xmax=2459
#punchim_ymax=1314
#windim_xmin=100
#windim_ymin=100
#windim_xmax=2362
#windim_ymax=2426
#thrshim_min=0
#thrshim_max=50
modeim_bin_size=1
modeim_kernel_width=15
'''
    image_params = get_image_params(imageset)

    # Set the LUNUS image parameters
    for pidx in range(len(image_params)):
        deck_and_extras = deck+image_params[pidx]
        A.LunusSetparamsim(pidx,deck_and_extras)
    A.LunusModeim()

    # Get the processed image
    lunus_filtered_data = flex.double()
    assert len(data)==self.n_panels # Jungfrau
    for pidx in range(len(data)):
      aye_panel = A.get_image_double(pidx)
      assert aye_panel.focus() == (254, 254)
      lunus_filtered_data.extend( aye_panel.as_1d() )
    lunus_filtered_data.reshape(flex.grid((self.n_panels,254,254)))
    self.view["lunus_filtered_data"] = lunus_filtered_data.as_numpy_array()

  def get_image_res_data(self):
      exp_data = self.expt.imageset.get_raw_data(0) # why a different access pattern? Are the data different?
      assert len(exp_data) == self.n_panels # 256 Jungfrau panels
      self.view["exp_data"] = [exp_data[pid].as_numpy_array() for pid in range(len(exp_data))]

  def write_hdf5(self,filenm):
      from simtbx.nanoBragg import utils
      img_sh = self.view["lunus_filtered_data"].shape
      assert img_sh == (self.n_panels,254,254)
      num_output_images = len([1 for key in self.view if self.view[key] is not None])
      print("Saving HDF5 data of shape %s to file %s"%(img_sh,filenm))
      beam_dict = self.expt.beam.to_dict()
      det_dict = self.expt.detector.to_dict()
      try:
        beam_dict.pop("spectrum_energies")
        beam_dict.pop("spectrum_weights")
      except Exception: pass
      import numpy as np
      comp = dict(compression='lzf' )
      with utils.H5AttributeGeomWriter(filenm,
                                image_shape=img_sh, num_images=num_output_images,
                                detector=det_dict, beam=beam_dict, dtype=np.float32,
                                detector_and_beam_are_dicts=True, compression_args=comp) as writer:
        for key in self.view:
          if self.view[key] is not None:
            writer.add_image(self.view[key])

def generate_phil_scope():
  from iotbx.phil import parse
  master_phil="""
    output {
      output_dir = .
        .type = path
      label = "lunus"
        .type = str
        .help = the label is used to construct the output file names, hdf5 and mask
      index = 0
        .type = int
        .help = serial number for output file
      enable = True
        .type = bool
        .help = whether to write an output image file (hdf5)
    }
    trusted_mask = ""
      .type = path
      .help = False values indicate bad pixels in the pickled flex.bool 3D grid
    refl = ""
      .type = path
      .help = The strong spot refl file (group A) to define the regions of interest
    expt = ""
      .type = path
      .help = The dials experiment file containing imageset, detector and beam
  """
  return parse(master_phil)
phil_scope = generate_phil_scope()

def parse_input():
  # The script usage
  import libtbx.load_env # implicit import
  help_message = """1) Use LUNUS assimilation filter to model the background image
2) Use the strong spots to construct a strong-spot ROI mask
Results can be viewed with dials.image_viewer <token>_%%%.hdf5 mask=<token>_%%%.mask
"""
  usage = ""
  from dials.util.options import OptionParser
  # Create the parser
  parser = OptionParser(usage=usage, phil=phil_scope, epilog=help_message)

  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)
  return params,options

def run(params):
    basename = "%s_%05d."%(params.output.label, params.output.index)
    M = mask_manager.from_files(params.trusted_mask, params.refl, params.expt)
    M.get_trusted_and_refl_mask()
    M.get_lunus_repl()
    M.get_image_res_data()
    M.write_hdf5(os.path.join(params.output.output_dir,basename+"hdf5"))
    M.resultant_mask_to_file(os.path.join(params.output.output_dir,basename+"mask"))

if __name__ == "__main__":
  params,options = parse_input()
  run(params)
