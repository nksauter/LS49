from __future__ import division
from six.moves import range
from scitbx.array_family import flex
from scitbx.matrix import sqr
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
from cctbx import crystal

"""Dev path to get step 4
develop a switch between time consuming SASE calc, and quick & dirty mono calc DONE
Check whether I can decrease crystal diffraction by factor of 10 DONE
Confirm orientations are random DONE
Sort out the output file names (make quick switch for intermediate file output) DONE
Put a record of the orientation in the file header. IT's IN THE SCREEN OUTPUT.
Test OpenMP on linux. DONE
"""
"""
problems: pixel values -40.  Fixed by using adc_offset_adu = 0
spot intensities flat:  (consequence of Tophat)
what is the crystal shape? Tophat vs Gauss vs Square makes a huge difference
water seemingly only diminishes bragg spots out to 2.76 Angstroms (consequence of direct transform res limit)
but it gives grey field throughout
how to specify a 500 um PAD ( NOT DONE, DOESN'T work , return to issue later)
how to specity Rayonix
apply_noise() seems to have added the psf even when add_psf() is off! --confirmed: observed when detector_psf_fwhm_mm!=0
There seem to be hot pixels in the noisified image!--values of 65535 or thereabouts (Fix by using adc_offset=10)
Why is the 'water ring' at 5.7 Angstroms?--Because of wavelength bug; use workaround
WHy does it look polarized in image_003 but not in noiseimage_001?
The size (cells_abc) seems to affect the intensity the wrong way --decreases instead of increases, but increases contrast over noise
    ---confirmed, image 1 seems to be normalized somehow
The size (cells_abs) affects radial width correctly (good)
C-centered symmetry not supported (non-90-degree cell angles not set)
"""

"""Initial try to produce simulations relevant to LS49.
step3:
produce a single PAD image
Ferredoxin, use Fcalc from PDB entry 1M2A, wild type ferredoxin from Aquifex aeolicus, 1.5 Angstrom
No anomalous signal at first.
Model a square PAD detector centered on direct beam, 2K x 2K, 0.11 mm pixels
Distance 141.7 mm, 500 um thick silicon.
Mosaicity angle is 0.05 stddev
500 mosaic domains
Domain size is 3000Angstrom (but it doesn't figure in to calculation).
flux is average 10^12 photons per shot, but is scaled by per-shot intensity
Illuminated volume is 3x3circlex10um long.
Choose a standard orientation at first, then a random orientation.
For the first shot, choose 7150 as our target beam energy.
Break it up into 100 spectral channels
Altogether this is 500 x 100 = 50000 mono-perfect subimages.
For the spectral dispersion, use shot #1, run 209 LG36.
Use a helium atmosphere, unattenuated beam passes through.
Work out mosaic rotation ensemble.
"""

"""Prospective plan
Simulate a single std-setting PAD image, using absolute units. (step 3)
Simulate a whole dataset, 20000 images.  Solve by molecular replacement (step 4)
Simulate an anomalous dataset, 7150 eV, but using fixed f',f". 100,000 images.  Solve by SAD.
*** Plug in the f' & f" as a function of wavelength.  Rerun the simulation at 7120 eV.
  Sort images into normalized-intensity wavelength bins and calculate maps
*** Put a detector at far distance, and one nearby, process both simultaneously.
  How to modify nanoBragg so it takes an offset Jungfrau, CSPAD, or Rayonix
  Make sure we can index both at the same time.
  Observe the energy positioning of Bragg spots, as function of deltapsi.
*** simulate polarized spectroscopy
"""

"""dials integration
dxtbx.print_header step4_000001.img
dials.import step4_000001.img # checks import
dials.estimate_gain datablock.json
dials.find_spots datablock.json threshold.dispersion.gain=1.47 filter.min_spot_size=2
dials.image_viewer datablock.json strong.pickle
dials.stills_process step4_00000[0-2].img threshold.dispersion.gain=1.47 filter.min_spot_size=2 indexing.known_symmetry.unit_cell=67.200,59.800,47.200,90.00,110.30,90.00 indexing.known_symmetry.space_group=C2
dials.image_viewer idx-step4_000000_integrated_experiments.json idx-step4_000000_integrated.pickle
"""

def tst_one(image,spectra,crystal,random_orientation):
  iterator = spectra.generate_recast_renormalized_image(image=image,energy=7150.,total_flux=1e12)

  quick = False
  if quick: prefix_root="step4Y_%06d"
  else: prefix_root="step4B_%06d"

  file_prefix = prefix_root%image
  rand_ori = sqr(random_orientation)
  from LS49.sim.step4_pad import run_sim2smv
  run_sim2smv(prefix = file_prefix,crystal = crystal,spectra=iterator,rotation=rand_ori,quick=quick)

if __name__=="__main__":
  from LS49.spectra.generate_spectra import spectra_simulation
  from LS49.sim.step4_pad import microcrystal
  SS = spectra_simulation()
  C = microcrystal(Deff_A = 4000, length_um = 1., beam_diameter_um = 1.0) # assume smaller than 10 um crystals
  mt = flex.mersenne_twister(seed=0)
  random_orientations = []
  N_total = 20160
  N_stride = 70 # total number of jobs
  for iteration in range(N_total):
    random_orientations.append( mt.random_double_r3_rotation_matrix() )

  import sys
  job_no = int(sys.argv[1])
  for idx in range(job_no,N_total,N_stride):
    print "idx------------------->",idx
    tst_one(image=idx,spectra=SS,crystal=C,random_orientation=random_orientations[idx])
  print "OK"
"""bsub -n 8 -o job${JOB}.log OMP_NUM_THREADS=8 libtbx.python ~/proj-1217/modules/LS49/sim/step4batch_pad.py ${JOB}
for JOB in `seq 10 69`; do bsub -q psanaq -n 8 -o job${JOB}.log OMP_NUM_THREADS=8 libtbx.python ~/proj-1217/modules/LS49/sim/step4batch_pad.py ${JOB};done

for JOB in `seq 62 69`; do bsub -q psanaq -o job${JOB}.log bsub -R "affinity[core(8)]" libtbx.python ~/proj-1217/modules/LS49/sim/step4batch_pad.py ${JOB};done FAIL
for JOB in `seq 62 69`; do bsub -q psanaq -n 24 -o job${JOB}.log OMP_NUM_THREADS=24 libtbx.python ~/proj-1217/modules/LS49/sim/step4batch_pad.py ${JOB};done
for JOB in `seq -w 0 69`; do bsub -q psanaq -n 12 -o jobA${JOB}.log OMP_NUM_THREADS=12 libtbx.python ~/proj-1217/modules/LS49/sim/step4batch_pad.py ${JOB};done
"""

