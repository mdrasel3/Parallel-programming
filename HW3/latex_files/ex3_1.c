#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// custom function my_Bcast_Loop
int my_Bcast_Loop(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {    
    int numprocs, rank;

    MPI_Comm_size(comm, &numprocs);
    MPI_Comm_rank(comm, &rank);

    if(rank == root) {
        for(int i = 0; i < numprocs; i++) {
            if(i != root) {
                MPI_Send(buffer, count, datatype, i, 0, comm);
            }
        }
    } else {
        MPI_Recv(buffer, count, datatype, root, 0, comm, MPI_STATUS_IGNORE);
    }

    return MPI_SUCCESS;
}

// main function
int main(int argc, char* argv[]) {

    int root = 0;
    int numprocs, rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    // allocating memory for an array of 100000 doubles
    double* data = (double*)malloc(sizeof(double)*100000);
    double* data_copy = (double*)malloc(sizeof(double)*100000);
    
    if(rank == root) {
        for(int i = 0; i < 100000; i++) {
            data[i] = i;
            data_copy[i] = i;

        }
    }
    
    // measuring and printing the time taken by MPI_Bcast
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    MPI_Bcast(data, 100000, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    printf("MPI_Bcast time(seconds) for rank %2d: %f\n", end-start,rank);
    
    // measuring and printing the time taken by MPI_Bcast_Loop
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    my_Bcast_Loop(data_copy, 100000, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    printf("my_Bcast_Loop time(seconds) for rank %2d: %f\n", end-start, rank);

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