from __future__ import division
from scitbx.array_family import flex
from six.moves import range
class csv:
  def __init__(self):
    self.energy = flex.double()
    self.red_fp = flex.double()
    self.red_fdp = flex.double()
    self.ox_fp = flex.double()
    self.ox_fdp = flex.double()
    with open("data_sherrell/FP_FDP_data.csv","rU") as F:
      lines = F.readlines()
      for line in lines:
        tokens = [float(f) for f in line.strip().split(",")]
        print tokens
        self.energy.append(tokens[0])
        self.red_fp.append(tokens[1])
        self.red_fdp.append(tokens[2])
        self.ox_fp.append(tokens[3])
        self.ox_fdp.append(tokens[4])
  def plot_them(self,plt,energies):
    keys = [ round(e,0) for e in energies if e >= 7070. and e<7201.]
    print keys
    elist = list(self.energy)
    print elist
    lookup = [ elist.index(key) for key in keys]
    plt.plot(elist, [ -self.red_fp[idx] for idx in lookup ], "b.")
    plt.plot(elist, [ self.red_fdp[idx] for idx in lookup ], "b.")
    plt.plot(elist, [ -self.ox_fp[idx] for idx in lookup ], "r.")
    plt.plot(elist, [ self.ox_fdp[idx] for idx in lookup ], "r.")

class george_sherrell:
  def __init__(self,file):
    self.energy = flex.double()
    self.fp = flex.double()
    self.fdp = flex.double()
    with open(file,"rU") as F:
      lines = F.readlines()
      for line in lines:
        tokens = [float(f) for f in line.strip().split()]
        self.energy.append(tokens[0])
        self.fp.append(tokens[1])
        self.fdp.append(tokens[2])
  def fp_fdp_at_wavelength(self,angstroms):
    lookup_energy = round(12398.425/angstroms,0)
    lookup_idx = list(self.energy).index(lookup_energy)
    return self.fp[lookup_idx], self.fdp[lookup_idx]
  def plot_them(self,plt,f1,f2):
    plt.plot(self.energy, self.fp, f1)
    plt.plot(self.energy, self.fdp, f2)


if __name__=="__main__":

  from LS49.spectra.generate_spectra import spectra_simulation
  SS = spectra_simulation()
  iterator = SS.generate_recast_renormalized_image(image=0,energy=7150.,total_flux=1e12)
  T = iterator.next()
  energies = list(12398.425/T[0])
  lambd = list(T[0])

  energies = flex.double(range(7000,7300))
  lambd = 12398.425/energies

  CSV = csv()

  wavelength = lambd[50]
  from iotbx import pdb
  pdb_text = open("1m2a.pdb","r").read()
  pdb_inp = pdb.input(source_info=None,lines = pdb_text)
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure.show_summary(prefix="Input structure ")
  #
  # take a detour to insist on calculating anomalous contribution of every atom
  scatterers = xray_structure.scatterers()
  for isc, sc in enumerate(scatterers):
    from cctbx.eltbx import sasaki, henke
    #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
    expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
    sc.fp = expected_henke.fp()
    sc.fdp = expected_henke.fdp()

    if sc.element_symbol() == "Fe":
      ife = isc
      break

  fp = flex.double()
  fdp = flex.double()
  for wave, energy  in zip(lambd,energies):
    print energy,
    sc = scatterers[ife]
    expected_henke = henke.table(sc.element_symbol()).at_angstrom(wave)
    sc.fp = expected_henke.fp()
    sc.fdp = expected_henke.fdp()

    print sc.fp, sc.fdp
    fp.append(sc.fp)
    fdp.append(sc.fdp)
    #from IPython import embed; embed()
  from matplotlib import pyplot as plt
  plt.plot (energies,fp,"r-")
  plt.plot (energies,fdp,"b-")
  CSV.plot_them(plt,energies)
  GS = george_sherrell("data_sherrell/pf-rd-ox_fftkk.out")
  GS.plot_them(plt,f1="b*",f2="g*")
  GS = george_sherrell("data_sherrell/pf-rd-red_fftkk.out")
  GS.plot_them(plt,f1="r*",f2="c*")
  GS = george_sherrell("data_sherrell/Fe.dat")
  GS.plot_them(plt,f1="m-",f2="y-")

  plt.show()
