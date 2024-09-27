#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=64
#SBATCH --partition=compute2011
#SBATCH --exclusive
##SBATCH --output=output.txt

module load mpi/openmpi/4.1.0
mpirun ./ex4_1
