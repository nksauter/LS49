
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from dxtbx.model.crystal import CrystalFactory
from scitbx.matrix import sqr

import numpy as np


def main(shape=shapetype.Tophat, cuda=False, seed=None):

    SIM = nanoBragg(verbose=10, oversample=0)
   
    # Defaults 
    cr = {'__id__': 'crystal',
          'real_space_a': (200, 0, 0),
          'real_space_b': (0, 180, 0),
          'real_space_c': (0, 0, 150),
          'space_group_hall_symbol': '-P 4 2'}
    cryst = CrystalFactory.from_dict(cr)
    SIM.Amatrix = sqr(cryst.get_A()).transpose().elems
        
    SIM.detpixels_fastslow = (64,64)
    SIM.pixel_size_mm=0.1
    SIM.wavelength_A=1.2
    SIM.verbose=10
    SIM.flux=1e12
    SIM.mosaic_spread_deg=0.02
    SIM.mosaic_domains=10
    SIM.polarization=1
    SIM.distance_mm=100
    SIM.F000=3e3  # FIXME this has to be equivalent to default_F, or else set in Fhkl, otherwise the test fails
    SIM.default_F=3e3
    SIM.progress_meter=True
    SIM.beamsize_mm=0.005
    SIM.exposure_s=1
    SIM.Ncells_abc=(15,15,15)
    SIM.show_params()

    # variable 
    SIM.xtal_shape=shape
    SIM.show_params() 

    if seed is not None:
        SIM.seed = seed
        SIM.randomize_orientation()
    if cuda:
        SIM.add_nanoBragg_spots_cuda()
    else:
        SIM.add_nanoBragg_spots()
    img = SIM.raw_pixels.as_numpy_array()
    return img


if __name__ == "__main__":
    failures = 0
    shapes = shapetype.Tophat, shapetype.Gauss, shapetype.Square, shapetype.Round
    for shape in shapes:
        img_cuda = main(cuda=True)
        img = main(cuda=False)
        if not np.allclose(img, img_cuda, atol=0.25):
            failures += 1 
    print ("There are %d fails" % failures)
    assert(failures == 0)
    print("OK")

