The file summit_make_reference.sh produces the comparison set used to
validate the ADSE13-196 simulations.  Note the following options used:

The GPU kernel code uses an exponential filter of -35; the reference
set however, uses the CPU add_nanoBragg_spots with no filter

The exascale_api is not used for the reference, in fact the GPU is not used.

The add spots algorithm is NKS, which other tests have shown to produce the
same numbers as JH or cuda.

The add_background() algorithm uses the stable_sort option, which also produces
the same results as the jh (original algorithm on CPU).  Slightly different
results +/- 1 come from the cuda algorithm.

The add_noise() function is commented out. Noise is commented out because
the current algorithm cannot be ported to CUDA giving same results.
