#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define N 65536

void gram_schmidt(double local_vectors[][N], double local_basis[][N], int local_m, int local_n, int rank) {
    int l, k;
    double t, norm;

    for (l = 0; l < local_m; l++) {
        t = local_vectors[l][rank * local_n];

        for (k = 0; k < l; k++) {
            double dot_product = 0.0;
            MPI_Allreduce(&t, &dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            for (int i = rank * local_n; i < (rank + 1) * local_n; i++) {
                t -= local_basis[k][i] * dot_product;
            }
        }

        norm = sqrt(t * t);
        for (int i = rank * local_n; i < (rank + 1) * local_n; i++) {
            local_basis[l][i] = t / norm;
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nproc = size;
    int local_n = N / nproc;
    double all_vectors[10][N];
    double local_vectors[10][N];
    double local_basis[10][N];

    if (rank == 0) {
        for (int i = 0; i < 10; i++) {
            char filename[20];
            sprintf(filename, "/work/korzec/LAB2/ex5/x_%d.dat", i);
            FILE *file = fopen(filename, "rb");
            fread(all_vectors[i], sizeof(double), N, file);
            fclose(file);
        }
    }

    MPI_Scatter(all_vectors, local_n * 10, MPI_DOUBLE, local_vectors, local_n * 10, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    gram_schmidt(local_vectors, local_basis, 10, N, rank);

    MPI_Gather(local_basis, N * 10, MPI_DOUBLE, all_vectors, N * 10, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double y[N];
        FILE *file_y = fopen("/work/korzec/LAB2/ex5/y.dat", "rb");
        fread(y, sizeof(double), N, file_y);
        fclose(file_y);

        double result[10] = {0.0};

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < N; j++) {
                result[i] += all_vectors[i][j] * y[j];
            }
        }

        printf("Dot products with y: ");
        for (int i = 0; i < 10; i++) {
            printf("%lf ", result[i]);
        }
        printf("\n");
    }

    MPI_Finalize();

    return 0;
}




