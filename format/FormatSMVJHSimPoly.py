"""Polychromatic SMV image reader for JHSim images."""
from __future__ import division

from dxtbx.format.FormatSMVJHSim import FormatSMVJHSim
from dxtbx.format.FormatSMV import FormatSMV
from dxtbx.model import Spectrum
from LS49.spectra.generate_spectra import spectrum_simulation

class FormatSMVJHSimPoly(FormatSMVJHSim):
  @staticmethod
  def understand(image_file):
    size, header = FormatSMV.get_smv_header(image_file)
    if header.get("BEAMLINE") != "fake": return False
    if header.get("PREFIX") is None: return False
    return header.get("PREFIX").find('batch')>=0

  def get_spectrum(self, index=0):
    ENERGY_CONV = 12398.419739640716
    img_prefix = self._header_dictionary["PREFIX"]
    idx = int(img_prefix[-6:])
    SS = spectrum_simulation()#common_data = common_data)
    SS.select(image=idx%100000)
    iterator = SS.generate_recast_renormalized_image(image=idx%100000,energy=7120.,total_flux=1e12)
    if False: SS.plot_recast_image(image=idx%100000, energy=7120.) # for debugging
    wavlen, flux, wavelength_A = next(iterator)
    energies = ENERGY_CONV/wavlen

    return Spectrum(energies, flux)

if __name__ == "__main__":
    import sys
    for arg in sys.argv[1:]:
        print(FormatSMVJHSimPoly.understand(arg))
