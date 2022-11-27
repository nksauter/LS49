from __future__ import division
from mpi4py import MPI
import os
mpi_rank = MPI.COMM_WORLD.Get_rank()
mpi_size = MPI.COMM_WORLD.Get_size()
print(mpi_rank, "of", mpi_size, "on node", os.environ.get("SLURMD_NODENAME",""))
