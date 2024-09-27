#include <mpi.h> 
#include <stdio.h>

void read_matrix(char* filename, double* A, int m, int N, int rank, int size, MPI_Comm comm);
void write_matrix(char* filename, double* A, int M, int N, int rank, int size, MPI_Comm comm);
void smear_matrix(double *A, double *B, int m, int N, double kappa, int rank, int size, MPI_Comm comm);

int main(int argc, char** argv) {

  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int M, N; // Global matrix dimensions
  int m = M/size; // Local block row size

  double *A, *B;

  // Read input matrix
  read_matrix("matrix.bin", A, m, N, rank, size, MPI_COMM_WORLD);

  // Perform smearing iterations
  for(int i=0; i<100; i++) {

    smear_matrix(A, B, m, N, 0.1, rank, size, MPI_COMM_WORLD);

    double *tmp = A;
    A = B; 
    B = tmp;

  }

  // Write output matrices
  char filename[100];
  sprintf(filename, "output_%d.bin", rank);
  write_matrix(filename, A, M, N, rank, size, MPI_COMM_WORLD);

  MPI_Finalize();
  return 0;

}

void read_matrix(char* filename, double* A, int m, int N, int rank, int size, MPI_Comm comm) {

  // Read local block row portion of matrix 
}

void write_matrix(char* filename, double* A, int M, int N, int rank, int size, MPI_Comm comm) {

  // Write local block row portion of matrix

}

void smear_matrix(double *A, double *B, int m, int N, double kappa, MPI_Comm comm) {

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Buffer to store rows being sent/received
  double *sendbuf, *recvbuf;

  // Allocate buffers
  sendbuf = malloc(sizeof(double) * N); 
  recvbuf = malloc(sizeof(double) * N);

  // Start non-blocking communication
  MPI_Irecv(recvbuf, N, MPI_DOUBLE, (rank-1+size)%size, 0, comm, &recv_req);
  MPI_Isend(sendbuf, N, MPI_DOUBLE, (rank+1)%size, 0, comm, &send_req);

  // Compute interior rows of B while communication happens
  for(int i=1; i<m-1; i++) {
    for(int j=0; j<N; j++) {
      B[i*N+j] = (1/(1+4*kappa))*(A[i*N+j] + kappa*(A[(i-1)*N+j] + A[(i+1)*N+j] + A[i*N+(j+1)%N] + A[i*N+(j-1+N)%N])); 
    }
  }

  // Wait for communication to finish
  MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
  MPI_Wait(&send_req, MPI_STATUS_IGNORE);

  // Copy boundary rows to buffers
  for(int j=0; j<N; j++) {
    sendbuf[j] = A[(m-1)*N+j]; 
    recvbuf[j] = A[0*N+j];
  }

  // Exchange boundary rows
  MPI_Sendrecv(sendbuf, N, MPI_DOUBLE, (rank+1)%size, 0, recvbuf, N, MPI_DOUBLE, (rank-1+size)%size, 0, comm, MPI_STATUS_IGNORE);

  // Compute boundary rows of B using received data
   for(int j=0; j<N; j++) {
    B[0*N+j] = (1/(1+4*kappa))*(recvbuf[j] + kappa*(A[1*N+j] + A[0*N+j+1] + A[0*N+j-1]));
    B[(m-1)*N+j] = (1/(1+4*kappa))*(sendbuf[j] + kappa*(A[(m-2)*N+j] + A[(m-1)*N+j+1] + A[(m-1)*N+j-1])); 
  }

  free(sendbuf);
  free(recvbuf);

}
