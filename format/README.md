## Additional dxtbx format classes for simulated images

### FormatSMVJHSimPoly.py
Read the LS49-style simulated images (fraom Sauter 2020 reference), but also make the spectrum available
to the programmer.  Installation and usage:
```
# make the format file available to dxtbx:
dxtbx.install_format -u FormatSMVJHSimPoly.py

# access the simulated LCLS spectrum within Python:
beam = experiment.beam
spectrum = experiment.imageset.get_spectrum(0)
energies_raw = spectrum.get_energies_eV().as_numpy_array()
weights_raw = spectrum.get_weights().as_numpy_array()
```
