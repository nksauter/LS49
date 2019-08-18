### Carryover results from CPU implementation, performed in April 2019

The benchmark runs the CCTBX/nanoBragg code to simulate
diffraction patterns from protein crystallography.  The output is 100,000 diffraction simulations from randomly oriented crystals. 
This represents 20 seconds data acquisition at a 5 kHz imaging rate, a capacity
expected at the LCLS-II-HE X-ray light source at SLAC in 2023.  The total size of the compressed images on disk is 1.00 TB (the images are
9 megapixels each).  

The calculation is embarassingly parallel: diffraction patterns are distributed over MPI 
worker ranks.  Within each diffraction image, the for loop over pixels is distributed 
over OpenMP threads.

Figure of Merit:  The total wall time for the simulation of 100,000 diffraction patterns, using 5% of the nodes on host computer.  This Table
will also be included in an upcoming:

| Host         | Elapsed |Nodes|MPI ranks |OpenMP threads/rank |Mean time/task|Mean time/rank=0 task|Effective throughput|
|--------------|---------|-----|----------|--------------------|--------------|---------------------|--------------------|
| edison       |  12.3 h | 280 |  6720    |  2                 | 2865 s       | 1950 s              | 2.26 Hz            |
| cori/haswell |  15.5 h | 120 |  3840    |  2                 | 1855 s       | 2062 s              | 1.79 Hz            |
| cori/knl     |   7.2 h | 484 |  8228    | 16                 | 1948 s       | 1991 s              | 3.86 Hz            |

### New results on Cori GPU node, August 2019

These prototype nodes each contain 8 Nvidia GPU cards.  There are 40 physical cores and 2 hyperthreads/core.  The results here
establish that we can effectively share 1 GPU card among 5 CPU cores each running an MPI worker rank. 
Load balancing works out for this particular application.
However, it is only feasible with this example to have 1 MPI rank / core; not 1 MPI rank / hyperthread

| Host         | Elapsed |Nodes|MPI ranks |Thread stride / rank |Kernel time|Inner loop time|Images out|Effective throughput|Scaled throughput|
|--------------|---------|-----|----------|-------------------|-----------|---------------|----------|--------------------|-----------------|
| cori/gpu     | 10 min  | 1   |   8      |  2                | 0.184 s   | 0.833 s       | 48       | 0.08 Hz            | 18.4 Hz         |
| cori/gpu     | 10 min  | 1   |   16     |  2                | 0.201 s   | 0.914 s       | 80       | 0.13 Hz            | 30.0 Hz         |
| cori/gpu     | 10 min  | 1   |   24     |  2                | 0.228 s   | 0.972 s       | 124      | 0.21 Hz            | 47.0 Hz         |
| cori/gpu     | 10 min  | 1   |   32     |  2                | 0.273 s   | 1.054 s       | 171      | 0.29 Hz            | 65.0 Hz         |
| cori/gpu     | 10 min  | 1   |   40     |  2                | 0.239 s   | 0.992 s       | 200      | 0.33 Hz            | 76.0 Hz         |
| cori/gpu     | 10 min  | 1   |   80     |  1                | 0.600 s   | 1.700 s       | 160      | 0.27 Hz            | 61.0 Hz         |

Explanation:  Each simulation is run on a single cori/gpu node.  Expanded resources are not available for the prototype machine.
However, there is every expectation that the application will weak-scale once the application is ported to Summit (image throughput
will increase proportionally to the # of nodes used).  Therefore, the "Scaled throughput" column estimates what the Summit throughput
will be if 5% of the GPU nodes are used (230 Summit nodes).

Detailed information:  within each image simulation the "inner loop" calculates diffraction contributions due to 100 separate X-ray 
energy channels.  The "inner loop time" listed in the table represents wall clock time for 1 of 100 iterations.  The "kernel time" 
is the subset of "inner loop time" that participates in host-to-device transfer, function evaluation, and device-to-host transfer.  


