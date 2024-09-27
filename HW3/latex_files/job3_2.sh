#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=48
#SBATCH --partition=compute2011
#SBATCH --exclusive
#SBATCH --output=output3_2.txt

module load mpi/openmpi/4.1.0
mpirun ./ex3_2

