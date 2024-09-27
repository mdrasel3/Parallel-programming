#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=64
#SBATCH --partition=compute2011
#SBATCH --exclusive

module load mpi/openmpi/4.1.0
mpirun ./ex5
