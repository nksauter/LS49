"""
Tests the pythony_beams setter in nanoBragg

The goal is to figure out how nanoBragg accumulates images over sources
and do it manually in exactly the same way... 

i.e. we want 1 N-beam call to nanoBragg to be the same as N 1-beam calls to nanoBragg

This test should be in simtbx once it passes.. if we stick with pythony_beams

"""
from simtbx.nanoBragg import shapetype
from dxtbx_model_ext import flex_Beam
from simtbx.nanoBragg import nanoBragg
from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from scitbx.matrix import sqr
import numpy as np

cr = {'__id__': 'crystal',
      'real_space_a': (300, 0, 0),
      'real_space_b': (0, 200, 0),
      'real_space_c': (0, 0, 150),
      'space_group_hall_symbol': '-P 4 2'}
cryst = CrystalFactory.from_dict(cr)

# nanoBragg property values as globals
AMAT = sqr(cryst.get_A()).transpose().elems
DEFAULT_F = 1e7
NCELLS_ABC = (15,15,15)
DET_SHAPE = (256,256)
POLA = 1
ROI = None
#ROI = ((90,92),(105,107))
#ROI = ((90,105),(105,120))


def test_xraybeams_laue(shape=shapetype.Square, insert_a_zero=False, cuda=False):
    """
    :param shape: nanoBragg crystal shape (nanoBragg.shapetype.Gauss, .Tophat, .Round, or .Square)
    :param insert_a_zero: whether to insert a 0 flux in one of the flex Beams
    :param cuda: whether to use cuda kernel
    :return:
    """
    Nbeams = 10  # number of beams (channels)
    banwd = 0.04  # xray bandwidth
    wavelens = np.linspace(
        1.2 - 1.2*banwd/2,
        1.2 + 1.2*banwd/2,
        Nbeams)
    np.random.seed(1)
    fluxes = np.random.uniform(0.5, 1.5, Nbeams) * 1e12
    
    # NOTE inserting a zero flux causes this test to fail
    if insert_a_zero:
        fluxes[0] = 0    

    # make a model beam
    beam_descr = {'direction': (-1, 0, 0),
                  'divergence': 0.0,
                  'flux': 1e12, 
                  'wavelength': 1e-10}   # overwrite wavelength and flux later

    xrbeams = flex_Beam()  # this can be appended to, right now its an empty beams list

    stepwise_spots = np.zeros(DET_SHAPE)
    for wav, fl in zip( wavelens, fluxes):

        nbr = nanoBragg(detpixels_slowfast=DET_SHAPE, verbose=10, oversample=0)
        if ROI is not None:
            nbr.region_of_interest = ROI
        nbr.default_F = DEFAULT_F
        nbr.xtal_shape = shape
        nbr.Ncells_abc = NCELLS_ABC
        nbr.wavelength_A = wav
        nbr.flux = fl
        nbr.polarization = POLA

        if cuda:
            nbr.add_nanoBragg_spots_cuda()
        else:
            nbr.add_nanoBragg_spots()
        
        stepwise_spots += nbr.raw_pixels.as_numpy_array()

        # keep track of beams for single call to nanoBragg using xray_beams
        beam = BeamFactory.from_dict(beam_descr)
        beam.set_wavelength(wav*1e-10)  # need to fix the necessity to do this..
        beam.set_flux(fl)
        beam.set_direction((-1, 0, 0))   # this is the convention, stick with it
        xrbeams.append(beam)

    nbr = nanoBragg(detpixels_slowfast=DET_SHAPE, verbose=10, oversample=0)
    if ROI is not None:
        nbr.region_of_interest = ROI
    nbr.xtal_shape = shape
    nbr.default_F = DEFAULT_F
    nbr.polarization = POLA
    nbr.Ncells_abc = NCELLS_ABC 
    nbr.xray_beams = xrbeams
    
    if cuda:  # TODO: test the cuda kernel sources mode, likely will need patches like the CPU one did
        nbr.add_nanoBragg_spots_cuda()
    else:    
        nbr.add_nanoBragg_spots()
    
    aggregate_spots = nbr.raw_pixels.as_numpy_array()

    assert(np.allclose(stepwise_spots, aggregate_spots))


if __name__ == "__main__":
    test_xraybeams_laue(shape=shapetype.Square)  # PASSES
    test_xraybeams_laue(shape=shapetype.Gauss)  # PASSES
    test_xraybeams_laue(shape=shapetype.Round)  # PASSES
    #test_xraybeams_laue(shape=shapetype.Tophat)  # FIXME: tophat is breaking!
    print("OK")

