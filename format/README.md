## Additional dxtbx format classes for simulated images

### FormatSMVJHSimPoly.py
Read the LS49-style simulated images (from Sauter 2020, Acta Crystallogr D76, 176-192), but also make the spectrum available
to the programmer.  Installation and usage:
```
# make the format file available to dxtbx:
dxtbx.install_format -u FormatSMVJHSimPoly.py

# access the simulated LCLS spectrum within Python:
from dxtbx.model.experiment_list import ExperimentListFactory
experiments = ExperimentListFactory.from_json_file(experiment_file, check_format=True)
experiment = experiments[0]
beam = experiment.beam
spectrum = experiment.imageset.get_spectrum(0)
energies_raw = spectrum.get_energies_eV().as_numpy_array()
weights_raw = spectrum.get_weights().as_numpy_array()
```
