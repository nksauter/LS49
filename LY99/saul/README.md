Planning for the March 2022 LY99 SPREAD data collection

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
