## SPREAD analysis of the March 2022 LY99 data.  

Development of an entirely new workflow, building on the Sauter 2020 and Mendez 2021 papers. 
<list of computational steps will be added here>
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

### 1. Optimize the detector metrology
### 2. Quality control on the unit cell
Over all runs, the means should vary less than the standard deviation within one run.
Prepare a single covariance model to be used throughout for unit cell acceptance/rejection.

### 3. Conventional merging
### 4. Refinement of the atomic model
Using the conventionally merged data, re-optimize the best available PDB model against the present unit cell.
Make sure to allow ADPs for the metal atoms, sulfur and higher (use PHENIX).  Also, since the SPREAD 
data are collected across a K-edge, the scattering factors for that particular metal (for example, Fe) 
should be refined.

### 5. Run the annulus worker
Print out statistics on spotfinder spots (indexed.refl, refined.expt).
Assess the useful resolution range for SPREAD.
Remove lattices that do not meet the filter conditions, generating a new set of expt + refl.
Run the [slurm script, annulus2.sh](./annulus2.sh). 

### 6. Run the trumpet plot
Create the trumpet plot.  Optionally perform DIALS refinement prior to plotting the statistics.  Filter 
 the outliers and remove them, thus generating a new set of expt + refl.  
 
Caution: savepng=True is appropriate if a subset of data are processed; but output is quite large if the 
 entire data are used, so in that case set savepng=False.  Sorry, ignore this advice (1/31/23) as there is 
 no correct way to disable savepng.  As currently implemented, savepng=False writes the plot to the 
 screen, which is counterproductive if executing within SLURM.
 
PNG files can be converted to an animated GIF if ImageMagick is present:
```
convert -delay 12 -loop 0 *.png shift2_annulus_hits.gif
```
Run the [slurm script, trial5_dials_refine.sh](./trial5_dials_refine.sh). 

### 7. Convert to integrated shoeboxes
Use the "specific" merge worker to replace the small shoeboxes of indexed strong spots,
 with the corresponding large shoeboxes from integration refls.  This is easily done interactively:

```
salloc -A m3890 -C cpu -t 60 -N 1 --qos=interactive
SIT_PSDM_DATA=/pscratch/sd/p/psdatmgr/data/pmscr srun -n 32 -c 4 cctbx.xfel.merge $MODULES/psii_spread/merging/application/specific/trial_switch.phil
```
The used phil file must be modified as follows:
```
# must reflect the filtered lattices from the previous step:
input.path=/pscratch/sd/n/nksauter/LY99/trumpet_plot/trial5/out

# must reflect the dials.stills_process output:
specific.integration.refl=/pscratch/sd/n/nksauter/LY99/dwpaley/2531058/%d/idx-data_%05d_integrated.refl
specific.integration.expt=/pscratch/sd/n/nksauter/LY99/dwpaley/2531058/%d/idx-data_%05d_integrated.expt
```
 
### 8. Grid search on mosaic parameters

For each combination of Ncells and eta, refine the orientation with diffBragg and determine ΔR, Correlation(ΔR,ΔΨ), and <σZ>
 on the specifc SPREAD annulus.  Pick the best global values for subsequent work, in this case Na=48, Nc=24, eta=0.0512 deg:
```
job      eta   Na  Nc   ΔR      CorrelΔR,ΔΨ  mean(σZ)
2959146  .0384 24   8 0.64px   -23.9%        1.3310   12 min 4 nodes, 6-7 balance
2959183  .0384 24  16 0.65px   -13.9%        1.2602   14 min 4 nodes, 6-7 balance
2959184  .0384 24  24 0.67px   -14.1%        1.2507   15 min 4 nodes, 6-7 balance
2959185  .0384 48   8 0.68px   -14.6%        1.3208   17 min 4 nodes, 6-7 balance
2959199  .0384 48  16 0.71px     5.7%        1.2549   18 min 4 nodes, 6-7 balance
2959200  .0384 48  24 0.79px    12.4%        1.3293   20 min 4 nodes, 6-7 balance
2959201  .0384 96   8 0.84px    -8.1%        1.3377   20 min 4 nodes, 6-7 balance
2959203  .0384 96  16 0.87px    21.3%        1.3833   20 min 4 nodes, 6-7 balance
2959204  .0384 96  24 0.95px    25.8%        1.3389   21 min 4 nodes, 6-7 balance
2959206  .0512 24   8 0.59px   -37.0%        1.3202   12 min 4 nodes, 6-7 balance
2959207  .0512 24  16 0.63px   -19.1%        1.2701   14 min 4 nodes, 6-7 balance
2959211  .0512 24  24 0.66px   -16.8%        1.2578   15 min 4 nodes, 6-7 balance
2959212  .0512 48   8 0.61px   -24.9%        1.2848   16 min 4 nodes, 6-7 balance
2959213  .0512 48  16 0.68px    -3.5%        1.2417   18 min 4 nodes, 6-7 balance
2959215  .0512 48  24 0.67px    -0.4%        1.2351   18 min 4 nodes, 6-7 balance <-- pick this
2959216  .0512 96   8 0.73px   -20.1%        1.3209   20 min 4 nodes, 6-7 balance
2959218  .0512 96  16 0.74px     9.0%        1.2637   19 min 4 nodes, 6-7 balance
2959219  .0512 96  24 0.81px    15.1%        1.2764   21 min 4 nodes, 6-7 balance
```
The SLURM file provided assumes that all experiments are taken as input, although a large subset could also suffice.
Source the command line script trial8_2923740_array.sh to submit the large array of jobs.  
Note that this script example assumes tetragonal or hexagonal, thus NcellsA = NcellsB.  It would have to be changed
e.g., for orthorhombic (photosystem II).  Caution should be taken to not use too many node hours.
                                                                                      
I have some concern over whether the code is correct, for I can no longer reproduce the small values of ΔR and mean(σZ).

### 9. Exa3A.  Run the exa3A script again with energy channel count narrowed from 2048 to 256.

Using diffBragg first derivatives to refine crystal rotation; also refine background and G-scale.

Run [slurm script, 3724210.sh](./3724210.sh). Resulting model (.expt) and data (.refl) are used for the next step.

### 10. Prepare model structure factors.

Cache the base structure factors to a pickle file to speed up the run time 
of the S1 final step.  This [slurm script, 3762323.sh](./3762323.sh) will compute reference structure factors, write them to 
pickle, and then go ahead with the "S1" final step listed below.  The rationale is that the reference structure factors are
somewhat time intensive.  They can be computed once and re-used for repeat trials of S1.  The following phil parameters are
specific to an instance of writing structure factors:
```
exafel.static_fcalcs.path=$WORK/reference/mmo_static_fcalcs.pickle
exafel.static_fcalcs.whole_path=$WORK/reference/mmo_miller_array.pickle
exafel.static_fcalcs.action=write
exafel.debug.format_offset=0
exafel.debug.energy_offset_eV=0
exafel.debug.energy_stride_eV=4.00 # new structure factors must be calculated if the stride changes
```
The following are calculated:
```
...path = *_static_fcalcs.pickle: 
 An extremely large file (e.g. 600 MB) with two items 
 1. A flex.miller_index array covering P1 in the chosen resolution range
 2. A flex.complex_double array of structure factors keyed to the Miller indices H.
    The columns are F_H(energy, n_channels) for F(bulk) + F(non-metal atoms). Thus this is
    the energy-dependent portion of the calculation not dependent on the metal model.

...whole_path = *_miller_array.pickle:
 A large file (e.g. 25 MB) with one item:
 A cctbx.miller_array with complex structure factors in P1 (Fmodel) calculated at the lowest energy channel.
```
### 11. Final step (S1) of scattering factor refinement.

This [slurm script, 4709418.sh](./4709418.sh) will refine scattering factors for two
independent Fe sites on a 4 eV grid spanning the K-edge. 

This is a self-contained example.  Summary:  48 nodes, 192 ranks, 78 min. runtime.  <br>On 3044 lattices from shift 2, analyze
66069 shoeboxes from 4.0-2.5Å, containing 23,121,402 pixels.


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
 - This work is now complete.  Here is a working [Slurm script, 5928113.sh](./5928113.sh) to apply Kramers-Kronig to Step 11.
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
   - All this is now done. Specific behavior is localized in sw1.py, and cases are selected with phil parameter exafel.metal= choice. Here is a [Slurm script](./5928113.sh) showing the use of exafel.metal=MMO2.
     - MMO2: methane monooxygenase, treating the two Fe sites as chemically distinct (realistic)
     - MMO2: methane monooxygenase, treating the Fe sites as chemically identical (for debug)
     - PSII: not implemented
 - Currently the scattering factor refinement stops after 1 macrocycle.  Need to stabilize behavior and then extend to 3(?) macrocycles as in Sauter (2020) paper.
