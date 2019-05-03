import libtbx
import scitbx
import math
from simtbx.nanoBragg import nanoBragg
from scitbx.matrix import sqr, col
from scitbx.array_family import flex
from simtbx.nanoBragg import shapetype
from libtbx.development.timers import Profiler
from boost.python import streambuf # will deposit printout into dummy StringIO as side effect


class ChannelSimulator:
  def __init__(self, rotation,sfall_main, N,
               mosaic_domains=25,
               mosaic_spread_deg=0.05,
               SEED=1,
               randomize=False):
    """

    :param rotation:
    :param sfall_main:
    :param N:
    :param mosaic_domains:
    :param mosaic_spread_deg:
    :param SEED:
    :param randomize:
    """
    UMAT_nm = flex.mat3_double()
    mersenne_twister = flex.mersenne_twister(seed=0)
    scitbx.random.set_random_seed(1234)
    rand_norm = scitbx.random.normal_distribution(
      mean=0,
      sigma=mosaic_spread_deg * math.pi/180.)
    g = scitbx.random.variate(rand_norm)
    mosaic_rotation = g(mosaic_domains)
    for m in mosaic_rotation:
      site = col(mersenne_twister.random_double_point_on_sphere())
      UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )

    Amatrix_rot = (rotation * sqr(sfall_main.unit_cell().orthogonalization_matrix())).transpose()

    self.SIM = nanoBragg(
                detpixels_slowfast=(3000,3000),
                pixel_size_mm=0.11,
                Ncells_abc=(N, N, N),
                wavelength_A=1,  # default we will update later@
                verbose=0)

    self.SIM.seed = SEED
    self.SIM.adc_offset_adu = 10 # Do not offset by 40
    self.SIM.mosaic_domains = mosaic_domains  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
    self.SIM.mosaic_spread_deg = mosaic_spread_deg # interpreted by UMAT_nm as a half-width stddev
    self.SIM.distance_mm=141.7
    self.SIM.distance_mm=141.7
    self.SIM.set_mosaic_blocks(UMAT_nm)
    self.SIM.oversample=1
    self.SIM.polarization=1
    self.SIM.default_F=1e5  # TODO: in the future we will init the energy dependent F_HKL here
    #self.SIM.Fhkl = energy_independent_F
    self.SIM.Amatrix_RUB = Amatrix_rot
    self.SIM.xtal_shape=shapetype.Gauss # both crystal & RLP are Gaussian
    self.SIM.progress_meter=False
    self.SIM.exposure_s=1.0 # so total fluence is e12
    self.SIM.beamsize_mm=0.003 #cannot make this 3 microns; spots are too intense
    if randomize:
      self.SIM.random_orientation()
    temp=self.SIM.Ncells_abc
    print("Ncells_abc=",self.SIM.Ncells_abc)
    self.SIM.Ncells_abc=temp

    # FIXME: add the CUDA init script here
    #initialize_GPU_variables()
    self.raw_pixels = self.SIM.raw_pixels  # FIXME: this will be on GPU

  def add_channel_pixels(self, wavelength_A, flux, rank, algo="cuda"):
    """
    :param wavelength_A:
    :param flux:
    :param rank:
    :param algo:
    :return:
    """
    self.SIM.flux=flux
    self.SIM.wavelength_A = wavelength_A

    P = Profiler("nanoBragg C++ rank %d"%(rank))
    if algo == "NKS":
      self.SIM.add_nanoBragg_spots_nks(streambuf(StringIO()))
      self.raw_pixels = self.SIM.raw_pixels  # NOTE: this actually incremenents hack for now, because cuda doesnt add spots
    elif algo == "JH":
      self.SIM.add_nanoBragg_spots()
      self.raw_pixels = self.SIM.raw_pixels
    elif algo == "cuda":
      #self.SIM.add_nanoBragg_spots_cuda()
      # FIXME: currently does not accumulate on the GPU

      self.SIM.allocate_cuda()
      self.SIM.add_energy_channel_cuda()
      self.SIM.get_raw_pixels_cuda()
      self.raw_pixels += self.SIM.raw_pixels  # NOTE: will be on GPU
      self.SIM.deallocate_cuda()
    else:
      raise Exception("unknown spots algorithm '%s' " % algo)
