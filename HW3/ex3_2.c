#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

// custom function to broadcast data in a tree-like fashion
int my_Bcast_Tree(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int numprocs, rank;

    MPI_Comm_size(comm, &numprocs);
    MPI_Comm_rank(comm, &rank);

    // finding smallest integer k with numprocs <= 2^k
    int k = ceil(log2(numprocs));

    // loop from k down to 1
    for(int i = k; i >= 1; i--) {
        // all processes with rank mod 2^i = root: send data to process with rank rank+2^(i-1) (if it exists)
        if(rank % (int)pow(2, i) == root) {
            int dest = rank + pow(2, i-1);
            if(dest < numprocs) {
                MPI_Send(buffer, count, datatype, dest, 0, comm);
            }
        } 
        // processes with rank mod 2^(i-1) = root: receive data from process with rank rank-2^(i-1)
        else if(rank % (int)pow(2, i-1) == root) {
            MPI_Recv(buffer, count, datatype, rank - pow(2, i-1), 0, comm, MPI_STATUS_IGNORE);
        }
    }

    return MPI_SUCCESS;
}

// main function
int main(int argc, char* argv[]) {
    int root = 0;
    int numprocs, rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // allocating memory for an array of 100000 doubles
    double* data = (double*)malloc(sizeof(double)*100000);
    double* data_copy = (double*)malloc(sizeof(double)*100000);

    // if this is the root process, initialize the data array
    if(rank == root) {
        for(int i = 0; i < 100000; i++) {
            data[i] = i;
	    data_copy[i] = i;
        }
    }

    // measure and print the time taken by MPI_Bcast
    double start = MPI_Wtime();
    MPI_Bcast(data, 100000, MPI_DOUBLE, root, MPI_COMM_WORLD);
    double end = MPI_Wtime();
    printf("MPI_Bcast time(seconds) for rank %2d: %f\n", end-start,rank);

    // measure and print the time taken by my_Bcast_Tree
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    my_Bcast_Tree(data_copy, 100000, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    printf("my_Bcast_Tree time(seconds) for rank %2d:%f\n", end-start, rank);

    // for verification
    for (int i = 0; i < 100000; i++){
            if(data[i] != data_copy[i]) {
                    printf("Mismatch at index %2d: %f vs %f\n", i, data[i], data_copy[i]);
                    break;
            }
    }
    
    free(data);
    free(data_copy);

    MPI_Finalize();
    return 0;
}