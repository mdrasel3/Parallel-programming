#include <mpi.h>
#include <stdio.h>

#define MAX_DIRECTIONS 1000

int main(int argc, char* argv[]) {
    
    int rank, numprocs;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int x0 = rank % 8; // calculate x0
    int x1 = rank / 8; // calculate x1

    int neighbors[4]; // declare an array for neighbours
    neighbors[0] = (x0 + 1) % 8 + x1 * 8; // positive x0
    neighbors[1] = x0 + ((x1 + 1) % 8) * 8; // positive x1
    neighbors[2] = (x0 + 7) % 8 + x1 * 8; // negative x0
    neighbors[3] = x0 + ((x1 + 7) % 8) * 8; // negative x1

    int pos = 0;
    int value = 0;
    int direction; 

    int directions[MAX_DIRECTIONS]; // declare an array for directions
    int num_directions = 0;

    // if rank = 0, open the file, otherwise print error and abort
    if (rank == 0) {
        FILE* file = fopen("directions.dat", "r");
        if (file == NULL) {
            perror("Error opening file");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        while (fscanf(file, "%d", &directions[num_directions]) != EOF && num_directions < MAX_DIRECTIONS) {
            num_directions++;
        }

        fclose(file);

        if (ferror(file)) {
            perror("Error reading file");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    // broadcast the number of directions and directions to all processors
    MPI_Bcast(&num_directions, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(directions, num_directions, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < num_directions; i++) {
        direction = directions[i];
        int n = neighbors[direction];

        MPI_Bcast(&n, 1, MPI_INT, pos, MPI_COMM_WORLD);

        if (rank == pos) {
            value += rank;
            MPI_Send(&value, 1, MPI_INT, n, 0, MPI_COMM_WORLD);
        }

        if (rank == n) {
            MPI_Recv(&value, 1, MPI_INT, pos, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        pos = n; // update the position
    }

    if (rank == 0) {
        printf("Final value: %d\n", value); // print the final value
    }

    MPI_Finalize();
    return 0;
}