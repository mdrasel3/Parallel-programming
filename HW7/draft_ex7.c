#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void fox_algorithm(int n, int n_local, double *A_local, double *B_local, double *C_local, int my_rank, int q);

int main(int argc, char **argv) {
    int my_rank, num_procs;
    int n, n_local;
    double *A, *B, *C;
    double *A_local, *B_local, *C_local;
    FILE *fileA, *fileB, *fileC;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (my_rank == 0) {
        // Process 0 reads matrix A from file
        fileA = fopen("A.txt", "r");
        if (fileA == NULL) {
            fprintf(stderr, "Error opening file A.txt\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fscanf(fileA, "%d %d", &n, &n);
        n_local = n / (int)sqrt(num_procs);

        A = (double *)malloc(n * n * sizeof(double));
        for (int i = 0; i < n * n; i++) {
            fscanf(fileA, "%lf", &A[i]);
        }
        fclose(fileA);

        // Process 0 reads matrix B from file
        fileB = fopen("B.txt", "r");
        if (fileB == NULL) {
            fprintf(stderr, "Error opening file B.txt\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fscanf(fileB, "%d %d", &n, &n);
        B = (double *)malloc(n * n * sizeof(double));
        for (int i = 0; i < n * n; i++) {
            fscanf(fileB, "%lf", &B[i]);
        }
        fclose(fileB);

        // Process 0 allocates memory for local matrices C
        C = (double *)malloc(n * n * sizeof(double));
    }

    // Broadcast n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    n_local = n / (int)sqrt(num_procs);

    // Allocate memory for local matrices A, B, C
    A_local = (double *)malloc(n_local * n_local * sizeof(double));
    B_local = (double *)malloc(n_local * n_local * sizeof(double));
    C_local = (double *)malloc(n_local * n_local * sizeof(double));

    // Scatter matrix A to all processes
    MPI_Scatter(A, n_local * n_local, MPI_DOUBLE, A_local, n_local * n_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Scatter matrix B to all processes
    MPI_Scatter(B, n_local * n_local, MPI_DOUBLE, B_local, n_local * n_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Perform matrix multiplication using the Fox algorithm
    fox_algorithm(n, n_local, A_local, B_local, C_local, my_rank, (int)sqrt(num_procs));

    // Gather the result matrices C from all processes to process 0
    MPI_Gather(C_local, n_local * n_local, MPI_DOUBLE, C, n_local * n_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        // Process 0 writes matrix C to file
        fileC = fopen("C.txt", "w");
        if (fileC == NULL) {
            fprintf(stderr, "Error opening file C.txt\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fprintf(fileC, "%d %d\n", n, n);
        for (int i = 0; i < n * n; i++) {
            fprintf(fileC, "%lf\n", C[i]);
        }
        fclose(fileC);

        // Free allocated memory
        free(A);
        free(B);
        free(C);
    }

    free(A_local);
    free(B_local);
    free(C_local);

    MPI_Finalize();

    return 0;
}

void fox_algorithm(int n, int n_local, double *A_local, double *B_local, double *C_local, int my_rank, int q) {
    int sqrt_q = (int)sqrt(q);
    int row_rank = my_rank / sqrt_q;
    int col_rank = my_rank % sqrt_q;

    MPI_Comm row_comm, col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, row_rank, col_rank, &row_comm);
    MPI_Comm_split(MPI_COMM_WORLD, col_rank, row_rank, &col_comm);

    double *temp_A = (double *)malloc(n_local * n_local * sizeof(double));
    double *temp_B = (double *)malloc(n_local * n_local * sizeof(double));

    for (int step = 0; step < sqrt_q; step++) {
        // Broadcast A_local in the row
        MPI_Bcast(A_local, n_local * n_local, MPI_DOUBLE, step % sqrt_q, row_comm);

        // Perform local matrix multiplication
        for (int i = 0; i < n_local; i++) {
            for (int j = 0; j < n_local; j++) {
                for (int k = 0; k < n_local; k++) {
                    C_local[i * n_local + j] += A_local[i * n_local + k] * B_local[k * n_local + j];
                }
            }
        }

        // Shift B_local in the column
        int dest = (col_rank - row_rank + sqrt_q) % sqrt_q;
        int source = (col_rank + row_rank) % sqrt_q;
        MPI_Sendrecv_replace(B_local, n_local * n_local, MPI_DOUBLE, dest, 0, source, 0, col_comm, MPI_STATUS_IGNORE);
    }

    free(temp_A);
    free(temp_B);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
}
