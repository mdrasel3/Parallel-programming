#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
    int rank, numprocs;
    double start_time, end_time; // Variables to store start and end times
	
    double local_sum;
    double local_sum_sq; // the sum of squares of the numbers
    double local_min;
    double local_max;
    int num_elements_per_file = 100;
    int num_files = 64;
    
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
   

    start_time = MPI_Wtime(); // Record the start time

    // Process data files based on rank
    for (int i = rank; i < num_files; i += numprocs) {
        char filename[20];
        snprintf(filename, sizeof(filename), "data_%d.dat", i);
        FILE *file = fopen(filename, "r");
        if (file == NULL) {
            fprintf(stderr, "Error: Could not open file %s\n", filename);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        // Process data from the file
        for (int j = 0; j < num_elements_per_file; ++j) {
            double num;
            if (fscanf(file, "%lf", &num) != 1) {
                fprintf(stderr, "Error: Invalid data format in file %s\n", filename);
                fclose(file);
                MPI_Abort(MPI_COMM_WORLD, -1);            
	    }
            // Update local calculations
            local_sum += num;
            local_sum_sq += num * num;
            if (num < local_min) local_min = num;
            if (num > local_max) local_max = num;
        }
        fclose(file);
    }

    // Perform reduction to calculate global statistics
    double global_sum, global_sum_sq, global_min, global_max;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Calculate mean and variance
    double mean = global_sum / (num_elements_per_file * num_files);
    double variance = (global_sum_sq / (num_elements_per_file * num_files)) - (mean * mean);
    // Print results
    if (rank == 0) {
        printf("Mean: %.2lf\n", mean);
        printf("Variance: %.2lf\n", variance);
        printf("Minimum: %.2lf\n", global_min);
        printf("Maximum: %.2lf\n", global_max);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime(); // Record the end time
    // Print the execution time on the root process (rank 0)
    if (rank == 0) {
    	double elapsed_time = end_time - start_time;
	printf("Execution Time: %.6f seconds\n",elapsed_time);
    }

    MPI_Finalize();
    return 0;
}
