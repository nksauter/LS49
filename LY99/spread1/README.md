## SPREAD analysis of the March 2022 LY99 data.  

Development of an entirely new workflow, building on the Sauter 2020 and Mendez 2021 papers. 
<list of computational steps will be added here>

### 1. Prepare model structure factors.

Cache the base structure factors to a pickle file to speed up the run time 
of the S1 final step.  Script will be checked in soon.

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

#### Known issues
 - An earlier execution of the same exact script finished in 40 minutes (on an earlier alcc-recipes build).
 
#### Scientific next steps
 - Add an additional 3000 lattices from shift 3.

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
 
