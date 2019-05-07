"""
Tests the pythony_beams setter in nanoBragg
Its really weird, things get out of control when DEFAULT_F becomes large

This test should be in simtbx once it passes.. if we stick with pythony_beams

"""
from dxtbx_model_ext import flex_Beam
from simtbx.nanoBragg import nanoBragg
from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from scitbx.matrix import sqr
import numpy as np

cr = {'__id__': 'crystal',
      'real_space_a': (200, 0, 0),
      'real_space_b': (0, 180, 0),
      'real_space_c': (0, 0, 150),
      'space_group_hall_symbol': '-P 4 2'}
cryst = CrystalFactory.from_dict(cr)
AMAT = sqr(cryst.get_A()).transpose().elems

DEFAULT_F = 1e7

def test_xraybeams_laue(insert_a_zero=False, cuda=False):
    Nbeams = 10
    det_shape = (256,256)
    wavelens = np.linspace(
        1.2 - 1.2*0.04,
        1.2 + 1.2*0.04,
        Nbeams)
    np.random.seed(1)
    fluxes = np.random.uniform(0.5,1.5,Nbeams ) * 1e12
    
    # NOTE inserting a zero flux causes this test to fail
    if insert_a_zero:
        fluxes[0] = 0    

    # make a model beam
    beam_descr = {'direction': (1, 0, 0),
                  'divergence': 0.0,
                  'flux': 1e12, 
                  'wavelength': 1e-10}   # overwrite wavelen and flux later

    xrbeams = flex_Beam()
    total_flux = np.sum( fluxes)
    
    stepwise_spots = np.zeros( det_shape)
    weights = []
    for wav, fl in zip( wavelens, fluxes):

        nbr = nanoBragg(detpixels_slowfast=det_shape,verbose=10)
        nbr.default_F = DEFAULT_F 
        nbr.Ncells_abc = (5,5,5)
        nbr.wavelength_A = wav
        nbr.flux = fl
        nbr.polarization = 1
        weight = Nbeams * fl / total_flux
        if cuda:
            nbr.add_nanoBragg_spots_cuda()
        else:
            nbr.add_nanoBragg_spots()
        stepwise_spots += weight * nbr.raw_pixels.as_numpy_array()
        weights.append( weight)    

        # keep track of beams for single call to nanoBragg using xray_beams 
        beam = BeamFactory.from_dict(beam_descr)
        beam.set_wavelength(wav*1e-10)
        beam.set_flux(fl)
        xrbeams.append(beam)
        
    nbr = nanoBragg(detpixels_slowfast=det_shape,verbose=10)
    nbr.default_F = DEFAULT_F
    nbr.polarization = 1
    nbr.Ncells_abc = (15,15,15)
    nbr.xray_beams = xrbeams
    
    if cuda:
        nbr.add_nanoBragg_spots_cuda()
    else:    
        nbr.add_nanoBragg_spots()
    
    aggregate_spots = nbr.raw_pixels.as_numpy_array()
    assert(np.allclose(stepwise_spots, aggregate_spots))

def test_xraybeams_2color():
    # simulate from two photon sources
    waveA = 1.5
    fluxA = 1e12
    waveB = 1.1
    fluxB = 1.85e11
   
    Ncells = (15,15,15)

    beam_descr = {'direction': (1, 0, 0),
                  'divergence': 0.0,
                   'polarization': 1., 
                  'flux': fluxA,
                  'wavelength': waveA*1e-10}

    # create two fresh beams
    beamA = BeamFactory.from_dict(beam_descr)
    beamB = BeamFactory.from_dict(beam_descr)

    # set corresponding wavelength and flux
    # FIXME, xray_beams expects beams with wavelength in meters
    beamA.set_wavelength(waveA*1e-10)
    beamA.set_flux(fluxA)
    beamB.set_wavelength(waveB*1e-10)
    beamB.set_flux(fluxB)

    xrbeams = flex_Beam()
    xrbeams.append(beamA)
    xrbeams.append(beamB)

    detshape=(128,128)
    
    ###
    # A sim
    nbrA = nanoBragg(detpixels_slowfast=detshape, verbose=10)
    nbrA.Amatrix = AMAT
    nbrA.polarization=1
    nbrA.default_F = DEFAULT_F 
    nbrA.Ncells_abc = Ncells 
    nbrA.wavelength_A = waveA
    nbrA.flux = fluxA
    nbrA.add_nanoBragg_spots()
    spotsA = nbrA.raw_pixels.as_numpy_array()

    ###
    # B sim
    nbrB = nanoBragg(detpixels_slowfast=detshape, verbose=10)
    nbrB.Amatrix = AMAT
    nbrB.polarization=1
    nbrB.default_F = DEFAULT_F 
    nbrB.Ncells_abc = Ncells 
    nbrB.wavelength_A = waveB
    nbrB.flux = fluxB
    nbrB.add_nanoBragg_spots()
    spotsB = nbrB.raw_pixels.as_numpy_array()

    ###
    # A+B sim
    nbr = nanoBragg(detpixels_slowfast=detshape, verbose=10)
    nbr.Amatrix = AMAT
    nbr.default_F = DEFAULT_F 
    nbr.Ncells_abc = Ncells
    nbr.xray_beams = xrbeams
    nbr.add_nanoBragg_spots()
    spotsAB = nbr.raw_pixels.as_numpy_array()
    
    # are the results of spotsA + spotsB the same ? 
    weightA = 2*fluxA / (fluxA+fluxB)
    weightB = 2*fluxB / (fluxA+fluxB)

    spotsAB2 = spotsA*weightA + spotsB*weightB
    
    assert(np.allclose(spotsAB, spotsAB2))

if __name__ == "__main__":
    #test_xraybeams_2color()
    test_xraybeams_laue()
    print("OK")

