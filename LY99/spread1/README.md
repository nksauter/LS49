## SPREAD analysis of the March 2022 LY99 data.  

Development of an entirely new workflow, building on the Sauter 2020 and Mendez 2021 papers. 
<list of computational steps will be added here>

### 1. Prepare model structure factors.

Cache the base structure factors to a pickle file to speed up the run time 
of the S1 final step.  This [slurm script, 3762323.sh](./3762323.sh) will compute reference structure factors, write them to 
pickle, and then go ahead with the "S1" final step listed below.  The rationale is that the reference structure factors are
very time intensive.  They can be computed once and re-used for repeat trials of S1.  The following phil parameters are
specific to an instance of writing structure factors:
```
exafel.static_fcalcs.path=$WORK/reference/mmo_static_fcalcs.pickle
exafel.static_fcalcs.whole_path=$WORK/reference/mmo_miller_array.pickle
exafel.static_fcalcs.action=write
exafel.debug.format_offset=0
exafel.debug.energy_offset_eV=0
exafel.debug.energy_stride_eV=4.00 # new structure factors must be calculated if the stride changes
```

### 2. S1. Final step of scattering factor refinement.

This [slurm script, 4709418.sh](./4709418.sh) will refine scattering factors for two
independent Fe sites on a 4 eV grid spanning the K-edge. 

This is a self-contained example.  Summary:  48 nodes, 192 ranks, 78 min. runtime.  <br>On 3044 lattices from shift 2, analyze
66069 shoeboxes from 4.0-2.5Å, containing 23,121,402 pixels.

#### Instructions for running the example on Perlmutter
 - Assure the following environment is sourced:
```
    export CFSSRC=<path-to-alcc-recipes-software-install>
    module purge
    module load PrgEnv-gnu cpe-cuda cudatoolkit
    module load evp-patch # known issue workaround
    source $CFSSRC/alcc-recipes/cctbx/activate.sh
    export MODULES=$CFSSRC/alcc-recipes/cctbx/modules
    export BUILD=$CFSSRC/alcc-recipes/cctbx/build
```
 - The following code diffs are required:
 ```
 --- a/simtbx/diffBragg/src/diffBragg.cpp
+++ b/simtbx/diffBragg/src/diffBragg.cpp
@@ -1815,6 +1815,7 @@ void diffBragg::add_diffBragg_spots(const af::shared<size_t>& panels_fasts_slows
 
     Npix_to_model = panels_fasts_slows.size()/3;
     SCITBX_ASSERT(Npix_to_model <= Npix_total);
+    raw_pixels_roi = af::flex_double(Npix_to_model); // NKS, only way to correctly size & zero array
     double * floatimage_roi = raw_pixels_roi.begin();
 
     diffBragg_rot_mats();
diff --git a/xfel/merging/command_line/merge.py b/xfel/merging/command_line/merge.py
index 3404842c5b..54678d568e 100644
--- a/xfel/merging/command_line/merge.py
+++ b/xfel/merging/command_line/merge.py
@@ -44,6 +44,7 @@ class Script(object):
   def __init__(self):
     self.mpi_helper = mpi_helper()
     self.mpi_logger = mpi_logger()
+    self.common_store = dict(foo="hello") # always volatile, no serialization, no particular dict keys guaranteed
 
   def __del__(self):
     self.mpi_helper.finalize()
@@ -163,6 +164,7 @@ class Script(object):
     # Perform phil validation up front
     for worker in workers:
       worker.validate()
+      worker.__dict__["common_store"] = self.common_store
     self.mpi_logger.log_step_time("CREATE_WORKERS", True)
 
     # Do the work
 ```

#### Known issues and bugs
 - An earlier execution of the same exact script finished in 40 minutes (on an earlier alcc-recipes build).
 - The current macrocycle implementation in sw1.py duplicates code; macrocycles 1-3 are unrolled.  Presumably this could be
a while loop.
 - There is a matplotlib warning (about memory consumption) in out/rank_0.err.  Code should be redesigned and warning
eliminated.
 - There is a silent Kokkos error that kills the job (used to be verbose in a previous build).  Presumably in kokkos_finalize.
Must be fixed eventually; error is encountered near the end, after results are output.

#### Scientific next steps & controls
 - Add an additional 3000 lattices from shift 3.
 - Refine the PDB model against the LY99 merged data (with PHENIX) so that the Fcalc's are exactly consistent with 
the unit cell used for SPREAD.
 - Said Phenix refinement should allow ADPs for all metals, and refine the Fe scattering factors.
 - Positive control: Now that we have the entire workflow, run it on the simulated data to see if we can refine scattering
factors.  Do this as a function of lattice count and energy granularity to see where the breaking point is.
 - Negative control: Miller index permutation
 - Negative control: Shuffle the energy channels
 - Sensible spot recruitment is an important next step
 - Implement the "last files" *.h5 image output, to verify that the final Bragg spot model with scattering factors lines
      up with experiment.

#### Implementation of new feature: Kramers-Kronig restraints
 - Summary: Provide a phil-based option to remove the Sauter (2020) restraints on 𝚫f′ and 𝚫f″, replacing them with 
 physics restraints based on the dispersion relations.
 - The existing phil controls for the 𝚫f restraints are defined here:  https://github.com/nksauter/LS49/blob/670579382699537ede522fe6fc857144167d2080/ML_push/phil.py#L68
 Turn the restraints off by setting fp.mean and fdp.mean to None.
 - As an example, the current function call to apply restraints is in the file https://gitlab.com/cctbx/psii_spread/-/blob/85d3c46220f62a0f7033ed8c3dbddecc7effad3c/merging/application/annulus/new_global_fdp_refinery.py. <br>
 The call to restrain 𝚫f′ is on line 286, and the call to restrain 𝚫f″ is on line 297.  The called function is here:
 https://gitlab.com/cctbx/psii_spread/-/blob/85d3c46220f62a0f7033ed8c3dbddecc7effad3c/merging/application/annulus/new_global_fdp_refinery.py#L308. <br>
 - Presumably a new function would be implemented to restrain 𝚫f′ and 𝚫f″ to each other.  The existing code serves as a working
 example explaining what the return values are (f = loss term; g1 = flex.double of the gradient of the loss term with respect
 to each scattering factor parameter), and how the return values are unpacked and fed into the minimizer.
 - A great approach would be to work out a unit test first, then actually implement the code.
   - The test would simulate f″ based on a very simple model of the K-edge.
   - Then use the dispersion relations to calculate f′.
   - Then sample both of these curves with Gaussian noise to simulate experimental measurement of the two curves.
   - Then develop a restraint model, and optimize the parameters.  Presumably use automatic differentiation for first-derivatives.
   - Compare the optimized model to the initial ground truth (and pass the test based on a tolerance). Show result in matplotlib.
 
#### Generalization of the code for photosystem II
 - As currently written the program will run out of memory due to the size of the structure factor table.  Total structure factor
 count scales linearly with a) unit cell volume, b) volume of the reciprocal space annulus requested, and c) the number of energy 
 channels analyzed for scattering factors. Also d) the number of metal classes analyzed (2 Fe for MMO, 4 Mn for PSII).
   - Breaking examples could be worked out with LY99/MMO to help redesign the code.
   - A leading candidate is to compute the wavelength-dependent 𝚫F directly in GPU memory (saving time) and possibly on the fly 
 (condensing the memory usage).  But all the details still need to be worked out. 
   - An initial PSII run might work if limited to low resolution (8Å) and very granular energy (8 eV).
 - Somehow we need to create a framework for parsing the PDB that is aware of the metal sites.  
   - For example, MMO Fe sites are 601 and 602.
   - We could analyze the data using 601 and 602 as different classes, or both in the same class.
   - For the four Mn atoms in PSII, we could treat the NCS monomers as being in the same class or a different class.
   - What code should be provided to parse these cases, how should we switch options using phil, and what hooks should be provided
 to the main code (which doesn't know the details of every structure)?  Presumably a Python class object whose data contains 
 information about the metals, and which is specialized for every use case.
   - The current code (labels=601, labels=602) needs to be generalized (labels=601,602).
   - The class would have to set it own preset_starting_model.
