#!/bin/bash
#SBATCH -J drag-force
#SBATCH -o out.%j
#SBATCH -n 256             # total number of mpi tasks requested
#SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t 12:00:00        # run time (hh:mm:ss)

# the default python on stampede has an old version of scipy... we
# need a newer one with minimize_scalar.  so before submitting job,
# run these commands:
#
# module load intel/14.0.1.106
# module load python/2.7.6

ulimit -c 0
set -x

ibrun python drag-force-mla.py
