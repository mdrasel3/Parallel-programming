#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=48
#SBATCH --partition=compute2011
#SBATCH --exclusive
##SBATCH --output=output_ex8.txt

module load mpi/openmpi/4.1.0
mpirun ./ex8
