#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define N 65536
#define M 10
#define FILE_PATH "/work/korzec/LAB2/ex5/x_%d.dat"

void read_vector(float *vec, int id) {
    char filename[50];
    sprintf(filename, FILE_PATH, id);
    FILE *file = fopen(filename, "r");
    fread(vec, sizeof(float), N, file);
    fclose(file);
}

float dot_product(float *a, float *b, int n_bar) {
    float sum = 0.0;
    for (int i = 0; i < n_bar; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

void normalize(float *a, int n_bar) {
    float norm = sqrt(dot_product(a, a, n_bar));
    for (int i = 0; i < n_bar; i++) {
        a[i] /= norm;
    }
}

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n_bar = N / size;
    float *x = (float *)malloc(n_bar * sizeof(float));
    float *b = (float *)malloc(n_bar * sizeof(float));

    for (int l = 0; l < M; l++) {
        if (rank == 0) {
            read_vector(x, l);
        }
        MPI_Bcast(x, n_bar, MPI_FLOAT, 0, MPI_COMM_WORLD);

        for (int k = 0; k < l; k++) {
            float dot = dot_product(b, x, n_bar);
            MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            for (int i = 0; i < n_bar; i++) {
                x[i] -= b[i] * dot;
            }
        }

        normalize(x, n_bar);
        memcpy(b, x, n_bar * sizeof(float));
    }

    free(x);
    free(b);
    MPI_Finalize();
    return 0;
}
