#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=compute2011
#SBATCH --exclusive

module load mpi/openmpi/4.1.0
mpirun ./ex9
