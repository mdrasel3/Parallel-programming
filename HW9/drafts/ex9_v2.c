#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void read_matrix(const char* filename, double** local_matrix, int* M, int* N, int* local_m, MPI_Comm comm) {
    MPI_File fh;
    MPI_Status status;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Open the file for reading
    MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    // Read M and N from the file
    if (rank == 0) {
        MPI_File_read(fh, M, 1, MPI_UNSIGNED, &status);
        MPI_File_read(fh, N, 1, MPI_UNSIGNED, &status);
        printf("Matrix dimensions read from file: M = %u, N = %u\n", *M, *N); // Print M and N
    }

    // Broadcast M and N to all processes
    MPI_Bcast(M, 1, MPI_UNSIGNED, 0, comm);
    MPI_Bcast(N, 1, MPI_UNSIGNED, 0, comm);

    // Calculate the size of the local block
    *local_m = *M / size;

    // Allocate memory for the local matrix
    *local_matrix = (double*)malloc((*local_m) * (*N) * sizeof(double));

    // Create a derived datatype for the local block
    MPI_Datatype file_type;
    int sizes[2] = {*M, *N};
    int subsizes[2] = {*local_m, *N};
    int starts[2] = {rank * (*local_m), 0};
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &file_type);
    MPI_Type_commit(&file_type);

    // Set the file view and read the local block
    MPI_File_set_view(fh, 2 * sizeof(int), MPI_DOUBLE, file_type, "native", MPI_INFO_NULL);
    MPI_File_read_all(fh, *local_matrix, (*local_m) * (*N), MPI_DOUBLE, &status);

    // Close the file and free the datatype
    MPI_File_close(&fh);
    MPI_Type_free(&file_type);
}

void write_matrix(const char* filename, double* local_matrix, int M, int N, int local_m, MPI_Comm comm) {
    MPI_File fh;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Create a derived datatype for the local block
    MPI_Datatype file_type;
    int sizes[2] = {M, N};
    int subsizes[2] = {local_m, N};
    int starts[2] = {rank * local_m, 0};
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &file_type);
    MPI_Type_commit(&file_type);

    // Open the file for writing
    MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // Set the file view and write the local block
    MPI_File_set_view(fh, 2 * sizeof(int), MPI_DOUBLE, file_type, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, local_matrix, local_m * N, MPI_DOUBLE, &status);

    // Close the file and free the datatype
    MPI_File_close(&fh);
    MPI_Type_free(&file_type);
}

void smear_matrix(double* local_A, double* local_B, int m, int N, double kappa, MPI_Comm comm) {
    int rank, size, up, down;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int ndims = 2; 
    int dims[2] = {m,N};
    int periods[2] = {1,0};
    int reorder = 0;
    MPI_Comm cart_comm;
    int cart_rank;
    int coords[ndims];

    // create communictor with cartesian topology
   MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm);

   // get rank with respect to new communicator
   MPI_Comm_rank(cart_comm, &cart_rank);

   // get position in cartesian process grid
   MPI_Cart_coords(cart_comm, cart_rank, ndims, coords);

    MPI_Cart_shift(comm, 0, 1, &up, &down);
    printf("cart_rank %2i, communication in positive 0 direction. source %2i, dest %2i\n",cart_rank, up, down);

    // Allocate memory for the send and receive buffers
    double* send_up = (double*)malloc(N * sizeof(double));
    double* send_down = (double*)malloc(N * sizeof(double));
    double* recv_up = (double*)malloc(N * sizeof(double));
    double* recv_down = (double*)malloc(N * sizeof(double));

    // Copy the first and last rows to the send buffers
    memcpy(send_up, local_A, N * sizeof(double));
    memcpy(send_down, local_A + (m - 1) * N, N * sizeof(double));

    // Start non-blocking sends and receives
    MPI_Request requests[4];
    
    MPI_Isend(send_up, N, MPI_DOUBLE, up, 0, comm, &requests[0]);
    MPI_Isend(send_down, N, MPI_DOUBLE, down, 1, comm, &requests[1]);
    MPI_Irecv(recv_up, N, MPI_DOUBLE, up, 1, comm, &requests[2]);
    MPI_Irecv(recv_down, N, MPI_DOUBLE, down, 0, comm, &requests[3]);

    // Compute the inner rows of B
    for (int i = 1; i < m - 1; i++) {
        for (int j = 0; j < N; j++) {
            int up_idx = (i - 1) * N + j;
            int down_idx = (i + 1) * N + j;
            int left_idx = i * N + (j - 1 + N) % N;
            int right_idx = i * N + (j + 1) % N;
            local_B[i * N + j] = (1.0 / (1.0 + 4.0 * kappa)) * (local_A[i * N + j] +
                kappa * (local_A[up_idx] + local_A[down_idx] +
                local_A[left_idx] + local_A[right_idx]));
        }
    }

    // Wait for non-blocking communication to complete
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

    // Compute the first and last rows of B
    for (int j = 0; j < N; j++) {
        int up_idx = (m - 1) * N + j;
        int down_idx = j;
        int left_idx = (j - 1 + N) % N;
        int right_idx = (j + 1) % N;
        local_B[j] = (1.0 / (1.0 + 4.0 * kappa)) * (local_A[j] +
            kappa * (recv_up[j] + local_A[up_idx] +
            local_A[left_idx] + local_A[right_idx]));
        local_B[(m - 1) * N + j] = (1.0 / (1.0 + 4.0 * kappa)) * (local_A[(m - 1) * N + j] +
            kappa * (local_A[down_idx] + recv_down[j] +
            local_A[(m - 2) * N + left_idx] + local_A[(m - 2) * N + right_idx]));
    }

    // Free allocated memory
    free(send_up);
    free(send_down);
    free(recv_up);
    free(recv_down);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Matrix dimensions and local block size
    int M, N, local_m;
    double* local_A;
    double* local_B;

    // Read the matrix from the file
    read_matrix("/work/korzec/LAB2/ex9/matrix.bin", &local_A, &M, &N, &local_m, comm);

    // Allocate memory for the local B matrix
    local_B = (double*)malloc(local_m * N * sizeof(double));

    // Smearing parameter
    double kappa = 0.1;

    // Perform smearing iterations
    int iterations[] = {5, 20, 100};
    for (int i = 0; i < sizeof(iterations) / sizeof(iterations[0]); i++) {
        int iter = iterations[i];
        for (int j = 0; j < iter; j++) {
            smear_matrix(local_A, local_B, local_m, N, kappa, comm);
            // Copy the smeared matrix back to A for the next iteration
            memcpy(local_A, local_B, local_m * N * sizeof(double));
        }
        // Write the smeared matrix to a binary file
        char filename[256];
        sprintf(filename, "smeared_matrix_%d.bin", iter);
        write_matrix(filename, local_B, M, N, local_m, comm);
    }

    // Free allocated memory
    free(local_A);
    free(local_B);

    MPI_Finalize();
    return 0;
}
