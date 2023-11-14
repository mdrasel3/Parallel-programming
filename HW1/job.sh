#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=compute2011
#SBATCH --exclusive


module load mpi/openmpi/4.1.0
mpirun ./ex1_2
