#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int M, N, m;
    double* local_data;
    double* recv_top_row;
    double* recv_bottom_row;
    double* local_B;

    MPI_File fh;
    MPI_Datatype filetype;
    MPI_Status status;

    MPI_File_open(MPI_COMM_WORLD, "/work/korzec/LAB2/ex9/matrix.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        MPI_File_read(fh, &M, 1, MPI_UNSIGNED, &status);
        MPI_File_read(fh, &N, 1, MPI_UNSIGNED, &status);
    }

    MPI_Bcast(&M, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    m = M / size;

    local_data = malloc(m * N * sizeof(double));
    recv_top_row = malloc(N * sizeof(double));
    recv_bottom_row = malloc(N * sizeof(double));
    local_B = malloc(m * N * sizeof(double));

    int sizes[2] = {M, N}, subsizes[2] = {m, N}, starts[2] = {rank * m, 0};
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_set_view(fh, 2 * sizeof(int), MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

    MPI_File_read_all(fh, local_data, m * N, MPI_DOUBLE, &status);

    MPI_File_close(&fh);

    MPI_Type_free(&filetype);

    double kappa = 0.1;
    int num_iterations[] = {5, 20, 100};

    for (int i = 0; i < sizeof(num_iterations) / sizeof(num_iterations[0]); i++) {
        for (int iter = 0; iter < num_iterations[i]; iter++) {
            MPI_Comm cart_comm;
            int dims[2] = {size, 1};
            int periods[2] = {0, 0};
            MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

            int coords[2], up_rank, down_rank;
            MPI_Cart_coords(cart_comm, rank, 2, coords);
            MPI_Cart_shift(cart_comm, 0, 1, &up_rank, &down_rank);

            // Start non-blocking sends and receives for top and bottom rows
            MPI_Request send_requests[2], recv_requests[2];
            MPI_Isend(local_data, N, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &send_requests[0]);
            MPI_Isend(local_data + (m - 1) * N, N, MPI_DOUBLE, (rank + 1) % size, 1, MPI_COMM_WORLD, &send_requests[1]);
            MPI_Irecv(recv_top_row, N, MPI_DOUBLE, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, &recv_requests[0]);
            MPI_Irecv(recv_bottom_row, N, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &recv_requests[1]);

            // Compute the inner rows of B while communication is ongoing
            for (int i = 1; i < m - 1; i++) {
                for (int j = 0; j < N; j++) {
                    double sum = local_data[i * N + j] + kappa * (
                        local_data[((i + 1) % m) * N + j] +
                        local_data[((i - 1 + m) % m) * N + j] +
                        local_data[i * N + (j + 1) % N] +
                        local_data[i * N + (j - 1 + N) % N]
                    );
                    local_B[i * N + j] = sum / (1 + 4 * kappa);
                }
            }

            // Wait for non-blocking communication to complete
            MPI_Waitall(2, send_requests, MPI_STATUSES_IGNORE);
            MPI_Waitall(2, recv_requests, MPI_STATUSES_IGNORE);

            // Compute the first and last rows of B after communication is finished
            for (int j = 0; j < N; j++) {
                double sum = local_data[j] + kappa * (
                    recv_top_row[j] +
                    local_data[(j + 1) % N] +
                    local_data[(j - 1 + N) % N]
                );
                local_B[j] = sum / (1 + 4 * kappa);
            }

            for (int j = 0; j < N; j++) {
                double sum = local_data[(m - 1) * N + j] + kappa * (
                    local_data[((m - 2 + m) % m) * N + j] +
                    local_data[(m - 1) * N + (j + 1) % N] +
                    local_data[(m - 1) * N + (j - 1 + N) % N] +
                    recv_bottom_row[j]
                );
                local_B[(m - 1) * N + j] = sum / (1 + 4 * kappa);
            }

            // Swap the local_data and local_B pointers for the next iteration
            double *temp = local_data;
            local_data = local_B;
            local_B = temp;

            MPI_Comm_free(&cart_comm);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        char filename[256];
        sprintf(filename, "smeared_matrix_%d_iterations.bin", num_iterations[i]);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

        MPI_File_set_view(fh, 2 * sizeof(int), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

        MPI_File_write_all(fh, local_data, m * N, MPI_DOUBLE, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);
    }

    free(local_data);
    free(recv_top_row);
    free(recv_bottom_row);
    free(local_B);

    MPI_Finalize();
    return 0;
}
