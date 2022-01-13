from __future__ import division
from LS49.spectra.generate_spectra import spectra_simulation

def test(idx):
  SS = spectra_simulation()
  iterator = SS.generate_recast_renormalized_image(image=idx%100000,energy=7120.,total_flux=1e12)
  wavlen, flux, wavelength_A = next(iterator) # list of lambdas, list of fluxes, average wavelength
  # the lambdas are equally spaced in energy, but expressed as wavelength in Angstroms
  # the fluxes are given in units of photons
  for iline in range(len(wavlen)):
    print ("lambda= %.4f Angstrom, flux=%.0f"%(wavlen[iline],flux[iline]))

if __name__=="__main__":
  idx = int(input("Enter the image number:"))
  print("Image %d"%idx)
  test(idx)

