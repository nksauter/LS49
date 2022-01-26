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
8. An attempt (failed) to refine 10 images with diffBragg stage 1 (as a call to hopper_utils.refine).  In the worker spread_roi.py, comment in the call to ds1.  Then use this input script: [roi_mini.sh](./roi_mini.sh).  Two complaints are a) segfault when refinement of mosaic rotation is commented in, b) divergence from unit cell starting model.  Advise studying the [diffBragg API FAQ](https://github.com/cctbx/cctbx_project/tree/master/simtbx/diffBragg#apifaq) in detail. 
9. An attempt (in progress) to perform the same refinement with the exascale_api.  In the worker spread_roi.py, comment in the call to exa1.  Then use this input script: [exa_mini.sh](./exa_mini.sh).

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
