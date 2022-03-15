Planning for the March 2022 LY99 SPREAD data collection.  Develop an entirely new workflow for SPREAD analysis building on the Sauter 2020 and Mendez 2021 papers.  Perform the following steps:
1. Simulate a 100,000-image dataset, on a Rayonix form factor, with and without a point spread function (PSF).  See [slurm script 918365](./sim_918365.sh).  I will develop the scattering factor analysis first without the PSF, and later consider the PSF.  Also, in order for shot-to-shot spectra to be associated with these simulated images, a special dxtbx.format file [must be installed](../../format). 
2. Indexing and integration with dials.stills_process, naive, with [index1_922530.sh](./index1_922530.sh) and [phil file index1.phil](./index1.phil).
3. Due to a half-pixel mismatch between simulated diffraction and metadata, the Level-0 detector position must be refined.  Steps for this are 
  a) ```dials.combine_experiments ../922530/idx-LY99_MPIbatch_000*.img_refined.expt ../922530/idx-LY99_MPIbatch_000*.img_indexed.refl reference_from_experiment.detector=0``` to combine just 30000 images, 
  b) then dials.refine [refine_level0_mcd.phil](./refine_level0_mcd.phil) combined.* output.experiments=refined_mcd.expt output.reflections=refined_mcd.refl
  c) comparison of the input and output with dials.show to determine the refined beam position.
4. Indexing and integration with dials.stills_process, using the updated detector position, with [index2_927185.sh](./index2_927185.sh) and [index2.phil](./index2.phil).
5. Unit cell analysis using a covariance model.  Output all unit cells (tdata file) with [tdata_928007.sh](./tdata_928007.sh).  Then run the covariance analysis with ```uc_metrics.dbscan file_name=928005/ly99sim_30000.tdata space_group=C12/m1 feature_vector=a,b,c eps=0.20 write_covariance=928005.tt metric=L2norm show_plot=True```, ultimately writing the covariance file covariance_ly99sim_30000.pickle.
6. A conventional merging run with Friedel mates separate, with [merge_928123.sh](./merge_928123.sh). 
7. A straight run with the [annulus worker](./annulus_929171.sh) to determine the integration shoeboxes and number of pixels available for SPREAD.
8. A series of regression tests for basic capabilities leading to global parameter refinement.
   - [test_product_0](./test_product_0_1285245.sh): The input worker is asked to keep_imagesets and read_image_headers.  Then simply test the ability of a downstream worker
    to read spectra and raw data arrays from the imageset data. The total count
    reported in the main log should agree with the number of lattices passing the unit cell filter.
    Currently FAILS on 100 nodes, presumably exceeding the on-node memory.  Logs only capture 90014 of 97878 images, 1659 of 3200 ranks.  No traceback.  We do not have to make this test work, see test 00.
   - [test_product_00](./test_product_00_1385292.sh): The input worker is asked to keep_imagesets but not read_image_headers.  Then test the ability of a downstream worker
    to read spectra and raw data arrays from a newly constructed imageset. The total count
    reported in the main log should agree with the number of lattices passing the unit cell filter.
    Currently SUCCEEDS on 100 nodes, 100K images, in under 3 minutes.
   - [test_product_1](./test_product_1_1385692.sh): The worker defines the set of shoeboxes for subsequent spread analysis, and creates a new bool array
    in the reflections table to flag these shoeboxes.  RMSD statistics are reported, both for the spread shoeboxes and for the whole set of reflections.
    This test succeeds in 2m:42s with 100 nodes/100K images.
   - [test_product_2](./test_product_2_1404150.sh): This worker demonstrates per-spot parameter refinement.
    It reads the spectrum and pixel data in the downstream worker (with per-image destruction).  It defines
    the set of shoeboxes for subsequent spread analysis.  It creates a model consisting of the DIALS
    shoebox plane fit (Rossmann) plus the whitelist monochromatic Bragg simulated in the Kokkos_gpu execution
    space.  Then it executes a per-spot planar fit to the rectangular shoebox border.  RMSD and Z-score are
    evaluated before and after.
    This test succeeds, but load balancing is an issue:

|    Workers dispatched             | 1 node/ 1 K images    |10 nodes/ 10 K images     | 100 nodes/ 100 K images  |
|-----------------------------------|-----------------------|--------------------------|--------------------------|
|using the merge.py balance worker  | | 4m:50s, 17-49 images/rank | 7m:30s, 14-52 images/rank|
|not using the balance worker       | 3m:10s, 30-32 images/rank | 3m:40s, 28-32 images/rank | 5m:10s, 26-32 images/rank|
    This test succeeds with 10 nodes/10K images (5 min.), and 100 nodes/100K images(7.5 min).
   - [test_product_3](./test_product_3.sh): This worker demonstrates per-lattice parameter refinement.
    Shoebox plane fit is optimized, but here it is done as a single LBFGS minimization for the whole
    lattice.  A generic parameter optimizer is prototyped.
    The test succeeds with 10 nodes/10K images (5 min.).
   - [test_product_3A](./test_product_3A.sh): First demonstration of a lattice-wide optimized
    parameter, in this case the overall scale factor G. Also, this implements a polychromatic model 
    instead of a monochromatic approximation, although we still use wavelength-independent structure factors.
    There is a slight radial r.m.s.d. improvement, polychromatic (0.79 px) vs. DIALS (0.99px).
    Due to perlmutter degradation, only tested with 1 node, 1K images (4 min).
   - [test_product_3B](./test_product_3B.sh): Same G scale-factor refinement, but use the diffBragg simulator.
    Simulated shoebox pixels are closely proportional to those of nanoBragg, with C.C.>99% for most images, C.C.>97% for all.
    Radial r.m.s.d. improvement, is similar to that with nanoBragg: polychromatic (0.81 px) vs. DIALS (0.99px).
    Due to perlmutter degradation, only tested with 1 node, 1K images (19 min). The increase in wall time is
    1) evaluation of simulated pixels each iteration, not done in 3A; 2) intrinsic wall time increase due to each-cycle
    round trip in diffBragg.
   - [test_product_4](./test_product_4.sh): Refinement of G scale-factor and crystal rotation.
    Radial r.m.s.d. improvement is now polychromatic (0.72 px) vs. DIALS (0.99px).  Wall time for 4 nodes, 1K images is 100 min
    with 32 ranks.
   - [DiffBragg stage 1, roi_mini.sh](./roi_mini.sh): Refine images with diffBragg stage 1 (as a call to hopper_utils.refine).
     As implemented, it refines the G-scale factor, orientation, unit cell, and Ncells_abc.  Mosaic spread eta is fixed at the
     ground truth value (if phil flag is set to refine eta, results diverge slightly). More information about diffBragg:
     [diffBragg API FAQ](https://github.com/cctbx/cctbx_project/tree/master/simtbx/diffBragg#apifaq).
     Radial r.m.s.d. improvement here is polychromatic (0.79 px) vs. DIALS (0.99px).  Wall time for 8 nodes, 1K images is 16 min
     with 64 ranks.

10. An attempt (in progress) to perform the same refinement with the exascale_api.  In the worker spread_roi.py, comment in the call to exa1.  Then use this input script: [exa_mini.sh](./exa_mini.sh).

Begining Dec. 18, 2021, perform a new simulation with the Rayonix form factor
```
        -N  -n  N_SIM wall weather
907749-- 1  4   320    276s  3.44s  new detector dimensions, new atten(small diffs), no noise or PSF
907902-- 1  4   320    170s  2.11s  new detector dimensions, atten=F, noise=F, PSF=F, CUDA
907910-- 1  4   320    137s  1.70s  new detector dimensions, atten=F, noise=F, PSF=F, KOKKOS
907915-- 1  4   320    209s  2.60s  new detector dimensions, atten=T, noise=F, PSF=F, KOKKOS
907924-- 1  4   320    723s  9.03s  new detector dimensions, atten=T, noise=T, PSF=F, KOKKOS
907920-- 1  4   320   1656s 20.66s  new detector dimensions, atten=T, noise=T, PSF=T, KOKKOS

910805-- 4 16   1280  1715s 21.26s  new detector dimensions, atten=T, noise=T, PSF=T, KOKKOS
910806--64 256  20480 1722s 21.26s  new detector dimensions, atten=T, noise=T, PSF=T, KOKKOS

Dec 21, the system has been rebooted without exclusive occupancy, 1 rank per GPU.
917146-- 4 32   1280   975s 22.26s  new detector dimensions, atten=T, noise=T, PSF=T, KOKKOS
917155-- 4 64   1280   773s 35.67s  new detector dimensions, atten=T, noise=T, PSF=T, KOKKOS
917161--64 1024 20480  792s 36.51s  new detector dimensions, atten=T, noise=T, PSF=T, KOKKOS
917162--64 1024 100000 3412s 33.72s  new detector dimensions, atten=T, noise=T, PSF=T, KOKKOS
917302--64 1024 100000 3356s 33.73s  new detector dimensions, atten=T, noise=T, PSF=F, KOKKOS

remove the nividia-smi[ Consider this script as the main simulation result! ]
918365--64 1024 100000 3334s 33.71s  new detector dimensions, atten=T, noise=T, PSF=F, KOKKOS

dials.stills_process:
922530  10 320  99987  999s 3.20s   integration of data without PSF (psff), run 917302
922532  10 320  99971  988s 3.16s   integration of data with the PSF (psft), run 917162

again after beam center refinement, dials.stills_process:
927185  10 320  99975  994s 3.20s   integration of data without PSF (psff), run 917302
927187  10 320  99975  981s 3.16s   integration of data with the PSF (psft), run 917162

tdata runs
928007: psff tdata output
928005: psft tdata output.  Use the covariance file for this one.

merging
928123: 32 2048 97673 132s: psff merging, Friedel mates separate.
928132: 32 2048 98586 125s: psft merging, Friedel mates separate.

try the annulus worker:
run 929158, 2.5-2.9 angstrom
run 929171, 2.1-2.5 angstrom

will need N=4 nodes for this problem size.  test on 1/10 data with 1 node
1161807: requesting diffBragg stage 1 eta_abc refinement leads to segmentation fault

work in progress on exascale API:
exa_mini.sh
```
