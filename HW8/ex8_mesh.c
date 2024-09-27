#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void read_matrix(char *filename, int n_local_rows, int n_local_cols, double *matrix_local, int my_rank, int q, int rows, int cols);
//void circular_shift(int direction, int *coords, int *new_coords, int q);
//void perform_shifts(int (*shifts)[2], int my_rank, int q, MPI_Comm cart_comm, int rows, int cols, int n_local_rows, int n_local_cols, double *matrix_local, double *temp_matrix);
void write_matrix(char *filename, int n_local_rows, int n_local_cols, double *matrix_local, int my_rank, int q, int rows, int cols);

int main(int argc, char **argv) {
    int my_rank, num_procs;
    MPI_Comm cart_comm;
    int dims[2] = {4, 12};
    int periods[2] = {1, 1};
    int coords[2];
    int shifts_m;
    int shifts_n;
    int cart_rank;
    File *shifts_file;
    double *matrix, *matrix_local, *temp_matrix;
    
    MPI_Comm row_comm;
    MPI_Comm col_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (num_procs != dims[0] * dims[1]) {
        fprintf(stderr, "Number of processes must be 4 * 12 = 48\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // create communictor with cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);

    // get rank with respect to new communicator
    MPI_Comm_rank(cart_comm, &cart_rank);

    // get position in cartesian process grid
    MPI_Cart_coords(cart_comm, my_rank, 2, coords);

    if (my_rank == 0) {
        // Process 0 reads shifts.dat file
        shifts_file = fopen("/work/korzec/LAB2/ex8/shifts.dat", "r");
        if (shifts_file == NULL) {
            fprintf(stderr, "Error opening shifts.dat\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //for (int i = 0; i < 1000; i++) {
            //fscanf(shifts_file, "%d %d", &shifts[i][0], &shifts[i][1]);
         fscanf(shifts_file, "%d", &shifts_m);
         fscanf(shifts_file, "%d", &shifts_n);

        fclose(shifts_file);
    }

    // Broadcast shifts to all processes
    MPI_Bcast(&shifts_m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&shifts_n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Read matrix from file
    int rows = 400, cols = 600;
    int n_local_rows = rows / dims[0];
    int n_local_cols = cols / dims[1];
    double*matrix_local = malloc(n_local_rows * n_local_cols * sizeof(double));
    double*temp_matrix = malloc(n_local_rows * n_local_cols * sizeof(double));

    read_matrix("/work/korzec/LAB2/ex8/notstirred.txt", n_local_rows, n_local_cols, matrix_local, cart_rank, dims[0], rows, cols);

    // Perform circular shifts based on shifts.dat
    //perform_shifts(shifts, my_rank, dims[0], cart_comm, rows, cols, n_local_rows, n_local_cols, matrix_local, temp_matrix);
    MPI_Cat_Sub(cart_comm, int free_cord[2]{0,1},&row_comm);
    MPI_Cat_Sub(cart_comm, int free_cord[2]{1,0},&col_comm);

    // Write the final matrix to file
    write_matrix("transformed_matrix.txt", n_local_rows, n_local_cols, matrix_local, my_rank, dims[0], rows, cols);

    free(matrix_local);
    free(temp_matrix);

    MPI_Finalize();

    return 0;
}

void read_matrix(char *filename, int n_local_rows, int n_local_cols, double *matrix_local, int my_rank, int q, int rows, int cols) {
    if (my_rank == 0) {
        FILE *file = fopen(filename, "r");
        if (file == NULL) {
            fprintf(stderr, "Error opening %s\n", filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int n = rows * cols;
        double matrix = malloc(n * sizeof(double));

        for (int i = 0; i < n; i++) {
            fscanf(file, "%lf", &matrix[i]);
        }

        fclose(file);

        MPI_Scatter(matrix, n_local_rows * n_local_cols, MPI_DOUBLE, matrix_local, n_local_rows * n_local_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        free(matrix);
    } else {
        MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, matrix_local, n_local_rows * n_local_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

/*void circular_shift(int direction, int *coords, int *new_coords, int q) {
    switch (direction) {
        case 0:
            new_coords[0] = (coords[0] + 1) % q;
            new_coords[1] = coords[1];
            break;
        case 1:
            new_coords[0] = (coords[0] + q - 1) % q;
            new_coords[1] = coords[1];
            break;
        case 2:
            new_coords[0] = coords[0];
            new_coords[1] = (coords[1] + 1) % 12;
            break;
        case 3:
            new_coords[0] = coords[0];
            new_coords[1] = (coords[1] + 11) % 12;
            break;
        default:
            break;
    }
}
*/
/*void perform_shifts(int (*shifts)[2], int my_rank, int q, MPI_Comm cart_comm, int rows, int cols, int n_local_rows, int n_local_cols, double *matrix_local, double *temp_matrix) {
    int coords[2], new_coords[2];
    MPI_Cart_coords(cart_comm, my_rank, 2, coords);

    for (int i = 0; i < 1000; i++) {
        circular_shift(shifts[i][0], coords, new_coords, q);

        int dest_rank;
        MPI_Cart_rank(cart_comm, new_coords, &dest_rank);

        MPI_Sendrecv(matrix_local, n_local_rows * n_local_cols, MPI_DOUBLE, dest_rank, 0, temp_matrix, n_local_rows * n_local_cols, MPI_DOUBLE, dest_rank, 0, cart_comm, MPI_STATUS_IGNORE);

        for (int j = 0; j < n_local_rows * n_local_cols; j++) {
            matrix_local[j] = temp_matrix[j];
        }
    }
}
*/
void write_matrix(char *filename, int n_local_rows, int n_local_cols, double *matrix_local, int my_rank, int q, int rows, int cols) {
    if (my_rank == 0) {
        FILE *file = fopen(filename, "w");
        if (file == NULL) {
            fprintf(stderr, "Error opening %s\n", filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int n = rows * cols;
        double *matrix = (double *)malloc(n * sizeof(double));

        MPI_Gather(matrix_local, n_local_rows * n_local_cols, MPI_DOUBLE, matrix, n_local_rows * n_local_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        fprintf(file, "%d %d\n", rows, cols);
        for (int i = 0; i < n; i++) {
            fprintf(file, "%lf\n", matrix[i]);
        }

        fclose(file);

        free(matrix);
    } else {
        MPI_Gather(matrix_local, n_local_rows * n_local_cols, MPI_DOUBLE, NULL, 0, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
    }
}
