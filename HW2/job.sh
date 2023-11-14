#!/bin/bash
#SBATCH --partition=compute2011
#SBATCH --exclusive            

module load mpi/openmpi/4.1.0

# Compile the MPI program
mpicc -o ex2 ex2.c -lm

# Run the MPI program with 1 core
echo "Running with 1 core"
mpirun -np 1 ./ex2

# Run the MPI program with 2 cores
echo "Running with 2 cores"
mpirun -np 2 ./ex2

# Run the MPI program with 4 cores
echo "Running with 4 cores"
mpirun -np 4 ./ex2

# Run the MPI program with 64 cores
echo "Running with 64 cores"
mpirun -np 64 ./ex2


