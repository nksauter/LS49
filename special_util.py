
import pylab as plt
import time

from dxtbx.model.experiment_list import ExperimentListFactory
import numpy as np
from simtbx.nanoBragg.utils import ENERGY_CONV
from scitbx.array_family import flex

FBG_VS_STOL = flex.vec2_double([
    (0, 2.57), (0.0365, 2.58), (0.07, 2.8), (0.12, 5), (0.162, 8), (0.18, 7.32), (0.2, 6.75),
    (0.216, 6.75), (0.236, 6.5), (0.28, 4.5), (0.3, 4.3), (0.345, 4.36), (0.436, 3.77), (0.5, 3.17)])

###
### These are methods checked into the diffBragg branch of cctbx
###

def spots_from_pandas(pandas_frame, Famp, oversample=0,
                      cuda=False, device_Id=0, time_panels=True,
                      njobs=1,
                      exascale=False):
    from joblib import Parallel, delayed
    from simtbx.nanoBragg.utils import flexBeam_sim_colors

    df = pandas_frame

    print("Loading experiment models")
    expt_name = df.exp_name.values[0]
    El = ExperimentListFactory.from_json_file(expt_name, check_format=False)
    expt = El[0]
    print("Done loading models!")
    assert len(df) == 1
    Ncells_abc = tuple(map(lambda x: int(round(x)), df.ncells.values[0]))
    spot_scale = df.spot_scales.values[0]
    beamsize_mm = df.beamsize_mm.values[0]
    total_flux = df.total_flux.values[0]
    #oversample = df.oversample.values[0]

    # get the optimized spectra
    if "spectrum_filename" in list(df):
        spectrum_file = df.spectrum_filename.values[0]
        pink_stride = df.spectrum_stride.values[0]
        fluxes, energies = load_spectra_file(spectrum_file, total_flux=total_flux,
                                             pinkstride=pink_stride)
    else:
        fluxes = np.array([total_flux])
        energies = np.array([ENERGY_CONV/expt.beam.get_wavelength()])
    lam0 = df.lam0.values[0]
    lam1 = df.lam1.values[0]
    if lam0 == -1:
        lam0 = 0
    if lam1 == -1:
        lam1 = 1
    wavelens = ENERGY_CONV / energies
    wavelens = lam0 + lam1*wavelens
    energies = ENERGY_CONV / wavelens

    crystal = expt.crystal
    crystal.set_A(df.Amats.values[0])

    panel_list = list(range(len(expt.detector)))
    pids_per_job = np.array_split(panel_list, njobs)

    sim_args = {"CRYSTAL":crystal, "DETECTOR": expt.detector, "BEAM": expt.beam, "Famp": Famp,
                "fluxes": fluxes, "energies": energies, "beamsize_mm":beamsize_mm,
                "Ncells_abc": Ncells_abc, "spot_scale_override": spot_scale,
                "cuda": cuda, "oversample": oversample,
                "time_panels":time_panels, "device_Id": device_Id}

    if exascale:
        sim_args["include_background"] = False
        sim_args["include_spots"] = True
        results = exascale_sim(**sim_args)

    else:
        def main(pids):
            sim_args["pids"] = pids
            results = flexBeam_sim_colors(**sim_args)
            return results
        results = Parallel(n_jobs=njobs)(delayed(main)(pids_per_job[jid]) for jid in range(njobs))
        results = [result for job_results in results for result in job_results]
        pids, imgs = zip(*results)
        order = np.argsort(pids)
        #assert np.all(order == np.arange(len(pids)))
        results = np.array([imgs[i] for i in order])

    return results


def load_spectra_file(spec_file, total_flux=1., pinkstride=1, as_spectrum=False):
    wavelengths, weights = np.loadtxt(spec_file, float, delimiter=',', skiprows=1).T
    if isinstance(wavelengths, float) and isinstance(weights, float):
        # the file had one entry:
        wavelengths = np.array([wavelengths])
        weights = np.array([weights])
    if pinkstride > len(wavelengths) or pinkstride == 0:
        raise ValueError("Incorrect value for pinkstride")
    wavelengths = wavelengths[::pinkstride]
    weights = weights[::pinkstride]
    energies = ENERGY_CONV/wavelengths
    FLUXES = weights / weights.sum() * total_flux
    if as_spectrum:
        return list(zip(list(wavelengths), list(FLUXES)))
    else:
        return FLUXES, energies


def save_spectra_file(spec_file, wavelengths, weights):
    data = np.array([wavelengths, weights])
    np.savetxt(spec_file, data.T, delimiter=',', header="wavelengths, weights")


def image_data_from_expt(expt, as_double=True):
    iset = expt.imageset
    if len(iset) == 0:
        raise ValueError("imageset should have 1 shot")
    if len(iset) > 1:
        raise ValueError("imageset should have only 1 shot. This expt has imageset with %d shots" % len(iset))
    try:
        flex_data = iset.get_raw_data(0)
    except Exception as err:
        assert str(
            type(err)) == "<class 'Boost.Python.ArgumentError'>", "something weird going on with imageset data"
        flex_data = iset.get_raw_data()
    if not isinstance(flex_data, tuple):
        flex_data = (flex_data,)
    img_data = np.array([data.as_numpy_array() for data in flex_data])
    if as_double:
        img_data = img_data.astype(np.float64)
    return img_data


def save_model_to_image(expt, model, output_img_file, save_experiment_data=False):
    npanels = len(expt.detector)
    n_img = 1
    if save_experiment_data:
        n_img = 2
    panelX, panelY = expt.detector[0].get_image_size()
    from simtbx.nanoBragg.utils import H5AttributeGeomWriter

    with H5AttributeGeomWriter(filename=output_img_file, image_shape=(npanels, panelY, panelX), num_images=n_img,
                               beam=expt.beam, detector=expt.detector) as writer:
        writer.add_image(model)
        if save_experiment_data:
            exp_data = image_data_from_expt(expt)
            writer.add_image(exp_data)
    print("Wrote model to image %s" % output_img_file)


def save_numpy_mask_as_flex(numpymask, outfile):
    import pickle
    from dials.array_family import flex
    flexmask = tuple((flex.bool(m) for m in numpymask))
    with open(outfile, "wb") as f:
        pickle.dump(flexmask, f)



def exascale_sim(
  DETECTOR, BEAM, CRYSTAL=None, Famp=None, energies=None, fluxes=None,
  background_wavelengths=None, background_wavelength_weights=None,
  background_total_flux=None, background_sample_thick_mm=0.5,
  density_gcm3=1, molecular_weight=18,
  cuda=False, oversample=0, Ncells_abc=(50, 50, 50),
  mos_dom=1, mos_spread=0, beamsize_mm=0.001,
  crystal_size_mm=0.01, device_Id=0,
  verbose=0, default_F=0, interpolate=0, profile="gauss",
  spot_scale_override=None, time_panels=False,
  add_water=False, add_air=False, water_path_mm=0.005, air_path_mm=0,
  adc_offset=0, readout_noise=3, psf_fwhm=0, gain=1, mosaicity_random_seeds=None, include_background=True,
    include_spots=True):

  tstart = time.time()
  from simtbx.nanoBragg.nanoBragg_beam import NBbeam
  from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
  from simtbx.nanoBragg.sim_data import SimData
  from simtbx.nanoBragg.utils import get_xray_beams
  from scipy import constants
  import numpy as np
  ENERGY_CONV = 10000000000.0 * constants.c * constants.h / constants.electron_volt

  nbBeam = NBbeam()
  nbBeam.size_mm = beamsize_mm
  nbBeam.unit_s0 = BEAM.get_unit_s0()

  if include_spots:
    assert fluxes is not None and energies is not None
    wavelengths = ENERGY_CONV / np.array(energies)
    nbBeam.spectrum = list(zip(wavelengths, fluxes))

  nbCrystal = None
  if include_spots:
    nbCrystal = NBcrystal()
    nbCrystal.dxtbx_crystal = CRYSTAL
    nbCrystal.Ncells_abc = Ncells_abc
    nbCrystal.symbol = CRYSTAL.get_space_group().info().type().lookup_symbol()
    nbCrystal.thick_mm = crystal_size_mm
    nbCrystal.xtal_shape = profile
    nbCrystal.n_mos_domains = mos_dom
    nbCrystal.mos_spread_deg = mos_spread

  S = SimData()
  S.detector = DETECTOR
  S.beam = nbBeam
  S.crystal = nbCrystal
  S.using_cuda = cuda
  S.using_omp = False
  S.add_air = add_air
  S.air_path_mm = air_path_mm
  S.add_water = add_water
  S.water_path_mm = water_path_mm
  S.readout_noise = readout_noise
  S.gain = gain
  S.psf_fwhm = psf_fwhm
  S.include_noise = False

  if mosaicity_random_seeds is not None:
    S.mosaic_seeds = mosaicity_random_seeds

  if Famp is not None:
      assert device_Id == Famp.get_deviceID()
  S.instantiate_nanoBragg(verbose=verbose, oversample=oversample, interpolate=interpolate,
    device_Id=device_Id, default_F=default_F, adc_offset=adc_offset)

  SIM = S.D  # the nanoBragg instance
  if spot_scale_override is not None:
    SIM.spot_scale = spot_scale_override
  if Famp is not None and Famp.get_nchannels() != 1:
      raise ValueError("Famp should be for 1 channel")

  from simtbx.gpu import exascale_api
  gpu_simulation = exascale_api(nanoBragg=SIM)
  gpu_simulation.allocate_cuda() # presumably done once for each image

  from simtbx.gpu import gpu_detector as gpud
  gpu_detector = gpud(deviceId=SIM.device_Id, detector=DETECTOR,
                      beam=BEAM)
  gpu_detector.each_image_allocate_cuda()

  if include_spots:
      gpu_simulation.add_energy_channel_from_gpu_amplitudes_cuda(
        0, Famp, gpu_detector)
      per_image_scale_factor = 1./len(energies)
      gpu_detector.scale_in_place_cuda(per_image_scale_factor)

  if include_background:
    SIM.beamsize_mm = beamsize_mm
    if background_wavelength_weights is None:
        background_wavelength_weights = [1]
    if background_wavelengths is None:
        background_wavelengths = [BEAM.get_wavelength()]
    wavelength_weights = np.array(background_wavelength_weights)
    weights = wavelength_weights / wavelength_weights.sum() * background_total_flux
    spectrum = list(zip(background_wavelengths, weights))
    xray_beams = get_xray_beams(spectrum, BEAM)
    SIM.xray_beams = xray_beams
    SIM.Fbg_vs_stol = FBG_VS_STOL
    SIM.flux=sum(weights)
    SIM.amorphous_sample_thick_mm = background_sample_thick_mm
    SIM.amorphous_density_gcm3 = density_gcm3
    SIM.amorphous_molecular_weight_Da = molecular_weight
    gpu_simulation.add_background_cuda(gpu_detector)

  packed_numpy = gpu_detector.get_raw_pixels_cuda().as_numpy_array()
  gpu_detector.each_image_free_cuda()
  print("done free")
  tdone = time.time()-tstart
  if time_panels:
      time_per_panel = tdone / len(DETECTOR)
      print("Took %f seconds to simulate %d panels (%.4f seconds/panel)" %(tdone, len(DETECTOR), time_per_panel))

  return packed_numpy


def extract_background(panel_imgs):
    from scipy.ndimage import filters
    background_imgs = []
    for i, img in enumerate(panel_imgs):
        smoother = filters.gaussian_filter(img, sigma=0.5)
        print("extracting bg %d / %d" %(i+1, len(panel_imgs)))
        background = filters.median_filter(smoother, size=(10,10))
        background_imgs.append(background)
    background_imgs = np.array(background_imgs)
    return background_imgs


def determine_bkgrnd_scale(data, model, stride=100):
    """
    labels background pixels and selects every `stride` pixels for fitting data to model
    data : image data (np array),
    model : background model (np array)  same shape as data
    """
    is_bkgrnd_pix = np.zeros(data.shape, np.bool)
    print("labeling background pixels")
    for i_panel, panel_img in enumerate(data):
        if i_panel % 50==0:
            print("labeling background pixels %d / %d" %(i_panel+1, len(data)))
        is_bkgrnd_pix[i_panel] = label_background_pixels(panel_img, thresh=2.5, iterations=2)

    data_pix = data[is_bkgrnd_pix].ravel()[::stride]  # every 100th pix
    model_pix = model[is_bkgrnd_pix].ravel()[::stride]

    print("Optimizing background model to data...")
    def func(x):
        return np.sum((data_pix-x[0]*model_pix)**2)

    from scipy.optimize import minimize
    opt = minimize(func, x0=np.array([1]), method="Nelder-Mead")
    if opt.success:
        print("Success, background scale is %f" % opt.x[0])
    else:
        print("Failure!")
    return opt


def label_background_pixels(roi_img, thresh=3.5, iterations=1):
    """
    iteratively determine background pixels in a subimg
    """
    img_shape = roi_img.shape
    img1 = roi_img.copy().ravel()   # 1-D version
    background_pixels = None
    while iterations > 0:
        if background_pixels is None:
            outliers = is_outlier(img1, thresh)
            m = np.median(img1[~outliers])
            out_and_high = np.logical_and(outliers, img1 > m)
            background_pixels = ~out_and_high
        else:
            where_bg = np.where(background_pixels)[0]
            outliers = is_outlier(img1[background_pixels], thresh)
            m = np.median(img1[background_pixels][~outliers])
            out_and_high = np.logical_and(outliers, img1[background_pixels] > m)
            background_pixels[where_bg[out_and_high]] = False
        iterations = iterations - 1

    return background_pixels.reshape(img_shape)


def is_outlier(points, thresh=3.5):
    """http://stackoverflow.com/a/22357811/2077270"""
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
