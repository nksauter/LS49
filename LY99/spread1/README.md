***SPREAD analysis of the March 2022 LY99 data.  Development of an entirely new workflow, building on the Sauter 2020 and Mendez 2021 papers. 
<list of computational steps will be added here>

**1. Prepare model structure factors.

Script will be checked in soon.  Cache the base structure factors to a pickle file to speed up the run time 
of the S1 final step.

**2. S1. Final step of scattering factor refinement.

This [slurm script, 4709418.sh](./4709418.sh) will refine scattering factors for two
independent Fe sites on a 4 eV grid spanning the K-edge. 

