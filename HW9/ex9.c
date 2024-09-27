#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define kappa 0.1

void read_matrix(const char *filename, int *M, int *N, double **local_A, int *local_rows, int rank, int size);
void write_matrix(const char *filename, double *local_A, int local_rows, int N, int rank, int size);
void smear_matrix(double *local_A, double *local_B, int local_rows, int N, int rank, int size, MPI_Comm comm);

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int M, N;
    double *local_A;
    int local_rows;

    // Load matrix from file and get dimensions
    read_matrix("/work/korzec/LAB2/ex9/matrix.bin", &M, &N, &local_A, &local_rows, rank, size);
     
    // Perform smearing iterations
    int iterations[] = {5, 15, 80};
    int iter = 0;
	for (int i = 0; i < 3; i++) {
        for (int iter = 0; iter < iterations[i]; iter++) {
            double *local_B = (double *)malloc(local_rows * N * sizeof(double));
            smear_matrix(local_A, local_B, local_rows, N, rank, size, MPI_COMM_WORLD);
            // Copy result back to A for the next iteration
            memcpy(local_A, local_B, local_rows * N * sizeof(double));
            free(local_B);
        }
        // Write the final smeared matrix to file
        iter += iterations[i];
        char outFile[50];
        sprintf(outFile, "output_%d.bin", iter);
        write_matrix(outFile, local_A, local_rows, N, rank, size);
    }
   
    free(local_A);
    
    MPI_Finalize();

    return 0;
}

void read_matrix(const char *filename, int *M, int *N, double **local_A, int *local_rows, int rank, int size) {
    MPI_File fh;
    MPI_Datatype subblock;;
    MPI_Status status;

    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    if (rank == 0) {
        // Read M and N from file by rank 0
        MPI_File_read(fh, M, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read(fh, N, 1, MPI_INT, MPI_STATUS_IGNORE);
    }

    // Broadcast M and N to all processes
    MPI_Bcast(M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate local rows for each process
    *local_rows = *M / size;

    // Allocate memory for the local matrix
    *local_A = (double *)malloc((*local_rows) * (*N) * sizeof(double));
    int sizes[2], lsizes[2], starts[2];

    // Only one dimension (2D matrix)
    sizes[0] = *M;
    sizes[1] = *N;
    lsizes[0] = *M / size;  // local_rows
    lsizes[1] = *N;
    starts[0] = rank * (*M) / size;
    starts[1] = 0;

    MPI_Type_create_subarray(2, sizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &subblock);
    MPI_Type_commit(&subblock);
    MPI_File_set_view(fh,  2 * sizeof(int), MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);
    MPI_File_read_all(fh, *local_A, (*local_rows) * (*N), MPI_DOUBLE, &status);

    MPI_File_close(&fh);
}

void write_matrix(const char *filename, double *local_A, int local_rows, int N, int rank, int size) {
    MPI_File fh;
    MPI_Datatype subblock;
    MPI_Status status;

    int sizes[2], lsizes[2], starts[2];
    int M = local_rows*size;
    // Only one dimension (2D matrix)
    sizes[0] = M;
    sizes[1] = N;
    lsizes[0] = local_rows;  // local_rows
    lsizes[1] = N;
    starts[0] = rank * local_rows;
    starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &subblock);
    MPI_Type_commit(&subblock);
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    MPI_File_write_at(fh, 0, &M, 1, MPI_INT, &status);
    MPI_File_write_at(fh, 1*sizeof(int), &N, 1, MPI_INT, &status);
    MPI_File_set_view(fh,  2 * sizeof(int), MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, local_A, local_rows*N, MPI_DOUBLE, &status); 
   

    MPI_File_close(&fh);
}

void smear_matrix(double *local_A, double *local_B, int local_rows, int N, int rank, int size, MPI_Comm comm) {
    MPI_Request reqs[4];
    int prev_rank = (rank - 1 + size) % size;
    int next_rank = (rank + 1) % size;

    double *temp_recv_top = (double *)malloc(N * sizeof(double));
    double *temp_recv_bottom = (double *)malloc(N * sizeof(double));

    // Start non-blocking send and receive for row exchange
    MPI_Isend(&local_A[0], N, MPI_DOUBLE, prev_rank, 0, comm, &reqs[3]);
    MPI_Irecv(temp_recv_bottom, N, MPI_DOUBLE, next_rank, 0, comm, &reqs[0]);
    MPI_Isend(&local_A[(local_rows - 1) * N], N, MPI_DOUBLE, next_rank, 1, comm, &reqs[2]);
    MPI_Irecv(temp_recv_top, N, MPI_DOUBLE, prev_rank, 1, comm, &reqs[1]);

    // Compute rows while communication is ongoing
    for (int i = 1; i < local_rows - 1; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i * N + j;
            local_B[idx] = (1.0 / (1.0 + 4.0 * kappa)) *(local_A[idx] + kappa * (local_A[(i - 1) * N + j] + local_A[(i + 1) * N + j] + local_A[i * N + (j + 1) % N] + local_A[i * N + (j - 1 + N) % N]));
        }
    }

    // Wait for the communication to finish
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    // Compute boundary rows i 0 and m -1
    for (int j = 0; j < N; j++) {
        int i = 0;
        int idx0 = i * N + j;

        local_B[idx0] = (1.0 / (1.0 + 4.0 * kappa)) * (local_A[idx0] + kappa * (temp_recv_top[j] + local_A[(i + 1) * N + j] + local_A[i * N + (j + 1) % N] + local_A[i * N + (j - 1 + N) % N]));
        i = local_rows -1;
        int idx1 = i * N + j;
        local_B[idx1] = (1.0 / (1.0 + 4.0 * kappa)) * (local_A[idx1] + kappa * (local_A[(i - 1) * N + j] + temp_recv_bottom[j] + local_A[i * N + (j + 1) % N] + local_A[i * N + (j - 1 + N) % N]));
    }

}

