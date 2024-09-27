#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <omp.h>

// Function to read a matrix from a file
void read_matrix(char *filename, double **matrix, int *rows, int *cols) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read matrix dimensions from the first row of the file
    fscanf(file, "%d %d", rows, cols);

    // Allocate memory for the matrix
    *matrix = (double *)malloc((*rows) * (*cols) * sizeof(double));

    // Read matrix elements from the file
    for (int i = 0; i < *rows * *cols; i++) {
        fscanf(file, "%lf", &(*matrix)[i]);
    }

    fclose(file);
}

// Function to write a matrix to a file
void write_matrix(const char *filename, double *matrix, int rows, int cols) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Write matrix dimensions to the first row of the file
    fprintf(file, "%d %d\n", rows, cols);

    // Write matrix elements to the file
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%.2f ", matrix[j + cols * i]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}


int main(int argc, char **argv) {
    
 MPI_Init(&argc, &argv);

 int rank, nproc;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &nproc);

 double *A, *B, *C;
 int m, n_A, n_B, k;

    // Process 0 reads matrices A and B
 if (rank == 0) {
    read_matrix("/work/korzec/LAB2/ex6/A.txt", &A, &m, &n_A);
    read_matrix("/work/korzec/LAB2/ex6/B.txt", &B, &n_A, &n_B);

        k = n_B; // Number of columns in matrix B
    }

    // Broadcast and scatter to all processes
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_A, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);


int trd,nthreads,i,j,k
double a[m][n_A];
double b[n_B][k];
double c[m][k];
double start,end,Time;
//Create parallel region

start = omp_get_wtime()
#pragma omp parallel shared(a,b,c) private(trd,i,j,l)
{
    trd = omp_get_thread_num();
    
    //Initialize matrix
    #pragma omp for
    for(i = 0; i < m; ++i)
        for(j = 0; j < n_A; ++j)
            a[i][j] = i + j;
    
    #pragma omp for
    for(i = 0; i < n_B; ++i)
        for(j = 0; j < k; ++j)
            b[i][j] = i * j;
            
    #pragma omp for
    for(i = 0; i < m; ++i)
        for(j = 0; j < k; ++j)
            c[i][j] = 0;
    
    #pragma omp for
    for(j = 0; i < n_A; ++j)
        for(l = 0; l < k; ++l)
            c[i][j] += a[i][l]*b[l][j];
    
  end = omp_get_wtime();
  Time = end - start;  
  MPI_Finalize(); 
  return 0;
}

