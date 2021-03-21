from __future__ import division
from dials.array_family import flex
import pickle
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.development.timers import Profiler

def top_75_iterator():
  import re
  from LS49.adse13_187.case_data import lookup_repo
  for ikey, key in enumerate(lookup_repo):
    assert ikey==key
    item = lookup_repo[key]
    match = re.search("shot([0-9]*)",item)
    event_idx = int(match.groups()[0])
    yield event_idx

class mask_manager:
  def __init__(self,trusted_mask,refl_table,expt):
    self.trusted_mask = trusted_mask
    self.refl_table = refl_table
    self.expt = expt

  def mask_introspection(self):
    for midx in range(256):
      print(midx, (self.trusted_mask[midx].count(True)), (self.trusted_mask[midx].count(False)))
    #conclusion:  mask True==good pixels; mask False==bad pixels
    # invert mask:  True==bad; False==good
    # on the invert mask, we can set the shoebox pixels to False

  def get_trusted_and_refl_mask(self):
    mask = self.trusted_mask
    Z = self.refl_table
    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    size = len(Z)

    self.shoebox_mask=[]
    for midx in range(256):
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
      for islow in range(max(0,bbox[2]-3), min(slow_size,bbox[3]+3)):
        for ifast in range(max(0,bbox[0]-3), min(fast_size,bbox[1]+3)):
          false_field[islow*slow_size + ifast]=True
          N_true += 1
          if mask[ipanel][islow*slow_size + ifast] is False: N_bad_pizel+=1
    self.resultant = []
    for midx in range(256):
      newmask = mask[midx].__and__(self.shoebox_mask[midx])
      self.resultant.append(newmask)
      if False and midx<4: print ("Panel",midx, "good px", mask[midx].count(True), "shoebox px",
             self.shoebox_mask[midx].count(True), "result", newmask.count(True))
    print (N_true,"pixels were visited in the %d shoeboxes (with borders)"%size)
    print (N_bad_pizel,"of these were bad pixels")

  def plot_pixel_histograms(self):
    exp_data = self.expt.imageset.get_raw_data(0) # experimental data
    Z = self.refl_table
    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    size = len(Z)
    all_pixels = flex.double()
    for sidx in range(size): #loop through the shoeboxes
      ipanel = P[sidx]
      false_field = self.shoebox_mask[ipanel]
      slow_size = false_field.focus()[0]
      fast_size = false_field.focus()[1]
      bbox = S[sidx].bbox
      islow_limits = (max(0,bbox[2]-3), min(slow_size,bbox[3]+3))
      ifast_limits = (max(0,bbox[0]-3), min(fast_size,bbox[1]+3))
      for islow in range(islow_limits[0], islow_limits[1]):
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          all_pixels.append(exp_data[ipanel][islow*slow_size + ifast])
    from matplotlib import pyplot as plt
    #plt.hist(all_pixels, len(all_pixels)//100, range=(-10,20000), facecolor="orange")
    plt.hist(all_pixels, len(all_pixels)//100, facecolor="orange")
    #plt.ylim((-10,500))
    plt.ylim((-10,100))
    plt.show()

  def resultant_mask_to_file(self,file_name):
    all_done = tuple(self.resultant)
    with open(file_name,"wb") as M:
      pickle.dump(all_done,M)

  @classmethod
  def from_files(cls, trusted_mask_file, refl_file, expt_file):
    with open(trusted_mask_file,"rb") as M:
      mask = pickle.load(M)
    refl_table = flex.reflection_table.from_file(refl_file)
    expts = ExperimentListFactory.from_json_file(expt_file,
                                              check_format=True)
    return cls(mask,refl_table,expts[0])

  def get_lunus_repl(self,filenm):
    P = Profiler("LUNUS")
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
    assert len(data)==256 # Jungfrau
    for pidx in range(len(data)):
      aye_panel = A.get_image_double(pidx)
      assert aye_panel.focus() == (254, 254)
      lunus_filtered_data.extend( aye_panel.as_1d() )
    lunus_filtered_data.reshape(flex.grid((256,254,254)))
    self.lunus_filtered_data = lunus_filtered_data.as_numpy_array()

    # then write the data
    from simtbx.nanoBragg import utils
    if True: # params.write_output:
      if True: # params.write_experimental_data:
        exp_data = self.expt.imageset.get_raw_data(0) # why is this a different access pattern? Are the data different?
        assert len(exp_data) == 256 # Jungfrau panels
      img_sh = self.lunus_filtered_data.shape
      assert img_sh == (256,254,254)
      num_output_images = 3 # 1 + int(params.write_experimental_data)
      print("Saving exascale output data of shape", img_sh)
      beam_dict = self.expt.beam.to_dict()
      det_dict = self.expt.detector.to_dict()
      try:
        beam_dict.pop("spectrum_energies")
        beam_dict.pop("spectrum_weights")
      except Exception: pass

      with utils.H5AttributeGeomWriter(filenm,
                                image_shape=img_sh, num_images=num_output_images,
                                detector=det_dict, beam=beam_dict,
                                detector_and_beam_are_dicts=True) as writer:
        writer.add_image(self.lunus_filtered_data)
        P = Profiler("shoebox planes")
        self.modify_shoeboxes()
        del P
        writer.add_image(self.lunus_filtered_data)

        if True: # params.write_experimental_data:
            exp_data = [exp_data[pid].as_numpy_array() for pid in range(len(exp_data))]
            writer.add_image(exp_data)
        print("Saved output to file %s" % (filenm))

  def modify_shoeboxes(self, verbose=False): # and printing the shoeboxes in verbose mode
    exp_data = self.expt.imageset.get_raw_data(0) # experimental data
    Z = self.refl_table
    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    size = len(Z)

    for sidx in range(size): #loop through the shoeboxes
      ipanel = P[sidx]
      false_field = self.shoebox_mask[ipanel]
      slow_size = false_field.focus()[0]
      fast_size = false_field.focus()[1]
      bbox = S[sidx].bbox
      islow_limits = (max(0,bbox[2]-3), min(slow_size,bbox[3]+3))
      ifast_limits = (max(0,bbox[0]-3), min(fast_size,bbox[1]+3))
      # print out the res-data
      for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = exp_data[ipanel][islow*slow_size + ifast]
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
      print()
      # print out the lunus-repl shoebox
      for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = self.lunus_filtered_data[ipanel,islow,ifast]
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
      print()
      # now create a 2nd-order fit to the data.  First implementation, no weighting.
      from LS49.adse13_187.adse13_221.smooth_fit import replacement_pixels
      FIT=replacement_pixels(self, ipanel, islow_limits, ifast_limits, shoebox=S[sidx])
      # print out the fit shoebox
      for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = FIT.model_T(islow,ifast)
          self.lunus_filtered_data[ipanel,islow,ifast]=value # reset lunus array
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
      print()
      # print out the data minus background-fit shoebox
      for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = exp_data[ipanel][islow*slow_size + ifast]-FIT.model_T(islow,ifast)
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
      print("---")
      #input()

def multiple_cases():
  def get_any_case():
    from LS49.adse13_187.adse13_221.ad_hoc_run795_lookup import conversion
    for idx,item in enumerate(top_75_iterator()):
      if idx>0: exit() # quick check the first one
      run_no = 795 # only look at run number 795
      class Empty: pass
      E = Empty()
      E.trusted_mask_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/trusted_Py3.mask"
      E.expt_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split_c/split_%04d.expt"%item
      E.out_file = "top_event_%04d.mask"%idx
      E.hdf5_file = "top_exa_%04d.hdf5"%idx
      E.refl_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_%04d.refl"%(
        conversion[run_no][item]
        )
      yield E
  for E in get_any_case():
    print(E.expt_file)
    print(E.refl_file)
    M = mask_manager.from_files(E.trusted_mask_file, E.refl_file, E.expt_file)
    M.get_trusted_and_refl_mask()
    #M.plot_pixel_histograms()
    # compute lunus-repl image.  write res-data and lunus-repl to an HDF5 file
    M.get_lunus_repl(E.hdf5_file)
    #M.modify_shoeboxes()
    M.resultant_mask_to_file(E.out_file)

def single_case():
  def get_case_1():
      class Empty: pass
      E = Empty()
      E.trusted_mask_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/trusted_Py3.mask"
      E.refl_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_0309.refl"
      E.expt_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split_c/split_0648.expt"
      E.out_file = "XXXevent_648.mask"
      E.hdf5_file = "exa_648.hdf5"
      return E
  E = get_case_1()
  M = mask_manager.from_files(E.trusted_mask_file, E.refl_file, E.expt_file)
  M.get_trusted_and_refl_mask()
  # compute lunus-repl image.  write res-data and lunus-repl to an HDF5 file
  M.get_lunus_repl(E.hdf5_file)
  M.resultant_mask_to_file(E.out_file)

if __name__ == "__main__":
  #single_case()
  multiple_cases()
