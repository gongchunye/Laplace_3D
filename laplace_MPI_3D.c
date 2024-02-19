/****************************************************************
 * Laplace MPI 3D
 * William Fung and Duy Hoang
 * modified by Chunye Gong, 2024.02.19
 *******************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include <string.h>

#define COLUMNS      100
#define LAYERS       100         // depth perception
#define ROWS_GLOBAL  100         // this is a "global" row count
#define NPES           4       // number of processors
#define ROWS (ROWS_GLOBAL/NPES)  // number of real local rows

// communication tags
#define DOWN     100
#define UP       101   

#define MAX_TEMP_ERROR 0.01

double Temperature[ROWS+2][COLUMNS+2][LAYERS+2];
double Temperature_last[ROWS+2][COLUMNS+2][LAYERS+2];
double Temperature_vis[LAYERS][COLUMNS][ROWS];

void initialize(int npes, int myid);
void track_progress(int iter, double dtt);


int main(int argc, char *argv[]) {

    int i, j, k;
    int max_iterations;
    int iteration=1;
    double dt;
    struct timeval start_time, stop_time, elapsed_time;

    int        npes;                // number of PEs
    int        myid;           // my PE number
    double     dt_global=100;       // delta t across all PEs
    MPI_Status status;              // status returned by MPI calls
    FILE *f_bov;
    FILE *f_double;

    // the usual MPI startup routines
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    // Initilize the PLANE type *****
    MPI_Datatype PLANE;
    MPI_Type_vector(COLUMNS+2, LAYERS+2, (COLUMNS+2), MPI_DOUBLE, &PLANE);
    MPI_Type_commit(&PLANE);

    // verify only NPES PEs are being used
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    
    if (npes != NPES) {
        if(myid==0) {
            printf("This code must be run with %d PEs\n", NPES);
        }
        MPI_Finalize();
        exit(1);
    };

    // PE 0 asks for input
    if (myid==0) {
      max_iterations= 4000; 
        printf("Maximum iterations [100-4000]?\n");
        //scanf("%d", &max_iterations);
    }
    // bcast max iterations to other PEs
    MPI_Bcast(&max_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myid==0) gettimeofday(&start_time,NULL);
    
    initialize(npes, myid);
    
    while ( dt_global > MAX_TEMP_ERROR && iteration <= max_iterations ) {
        // main calculation: average my six neighbors
        for(i = 1; i <= ROWS; i++) {
            for(j = 1; j <= COLUMNS; j++){
                for(k = 1; k <= LAYERS; k++) {
					//Jacobi
                    // Temperature[i][j][k] = (Temperature_last[i+1][j][k] + Temperature_last[i-1][j][k] +
                                            // Temperature_last[i][j+1][k] + Temperature_last[i][j-1][k] +
                                            // Temperature_last[i][j][k+1] + Temperature_last[i][j][k-1])/6;
					//Gauss-Seidel
                    Temperature[i][j][k] = (Temperature_last[i+1][j][k] + Temperature[i-1][j][k] +
                                            Temperature_last[i][j+1][k] + Temperature[i][j-1][k] +
                                            Temperature_last[i][j][k+1] + Temperature[i][j][k-1])/6;
                }
            }
        }
        
        // COMMUNICATION PHASE: send and receive ghost planes for next iteration
        if (myid != npes - 1) {
            MPI_Send(&(Temperature[ROWS][0][0]), 1, PLANE, myid + 1, DOWN, MPI_COMM_WORLD);
        }
        if (myid != 0) {
            MPI_Recv(&(Temperature_last[0][0][0]), 1, PLANE, myid - 1, DOWN, MPI_COMM_WORLD, &status);
        }
        if (myid != 0) {
            MPI_Send(&(Temperature[1][0][0]), 1, PLANE, myid - 1, UP, MPI_COMM_WORLD);
        }
        if (myid != npes - 1) {
            MPI_Recv(&(Temperature_last[ROWS + 1][0][0]), 1, PLANE, myid + 1, UP, MPI_COMM_WORLD, &status);
        }

        dt = 0.0;

        for(i = 1; i <= ROWS; i++){
            for(j = 1; j <= COLUMNS; j++) {
                for(k = 1; k <= LAYERS; k++) {
                    //dt = fmax( fabs(Temperature[i][j][k]-Temperature_last[i][j][k]), dt);
					if (fabs(Temperature[i][j][k]-Temperature_last[i][j][k]) > dt){
						dt = fabs(Temperature[i][j][k]-Temperature_last[i][j][k]);
					}
                    Temperature_last[i][j][k] = Temperature[i][j][k];
                }
            }
        }

        // find global dt
        MPI_Allreduce(&dt, &dt_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // periodically print test values - only for PE in lower corner
        if((iteration % 300) == 0) {
            if (myid == npes-1){
                track_progress(iteration, dt_global);
	        }
        }

	iteration++;
    }

    // Slightly more accurate timing and cleaner output 
    MPI_Barrier(MPI_COMM_WORLD);

    // PE 0 finish timing and output values
    if (myid==0){
        gettimeofday(&stop_time,NULL);
	timersub(&stop_time, &start_time, &elapsed_time);

	printf("\nMax error at iteration %d was %f\n", iteration-1, dt_global);
	printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    }

    ////////// final bov files /////////////
    char values[] = ".values";
    char pe_num[20];
    snprintf(pe_num, sizeof(pe_num), "%d", myid); //itoa
    strcat(pe_num, values); // pe_num = pe_num + values
    char filename2[40] = "bovf";
    strcat(filename2, pe_num);
	//printf("%s\n",filename2);

    f_double = fopen(filename2, "wb");
    for (i = 1; i <= ROWS; i++) {
        for (j = 1; j <= COLUMNS; j++) {
            for (k = 1; k <= LAYERS; k++) {
                Temperature_vis[k - 1][j - 1][i - 1] = Temperature[i][j][k];
            }
        }
    }
    fwrite((void *)Temperature_vis, sizeof(double), ROWS*COLUMNS*LAYERS, f_double);
    fclose(f_double);

    char bov[] = ".bov";
    char pe_num2[20];
    snprintf(pe_num2, sizeof(pe_num2), "%d", myid); //itoa
    strcat(pe_num2, bov);
    char filename[40] = "laplacef";
    strcat(filename, pe_num2);
    f_bov = fopen(filename, "wb");
    if (f_bov == NULL) {
        fprintf(stderr, "Could not write data file \n");
        return -1;
    }
    fprintf(f_bov, "TIME: 100\n");
    fprintf(f_bov, "DATA_FILE: bovf%d.values\n", myid);
    fprintf(f_bov, "DATA_SIZE: %d %d %d\n", ROWS, COLUMNS, LAYERS);
    fprintf(f_bov, "DATA_FORMAT: DOUBLE\n");
    fprintf(f_bov, "VARIABLE: U\n");
    fprintf(f_bov, "DATA_ENDIAN: LITTLE\n");
    fprintf(f_bov, "CENTERING: ZONAL\n");
    fprintf(f_bov, "BRICK_ORIGIN: 0. 0. 0.\n");
    fprintf(f_bov, "BRICK_SIZE: 1.0 1.0 1.0\n");
    fclose(f_bov);
    /////////////////////////////////
    MPI_Type_free(&PLANE);
    MPI_Finalize();
}

void initialize(int npes, int myid){

    //double tMin, tMax;  //Local boundary limits
    int i,j,k;
    
    for(i = 0; i <= ROWS+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            for(k = 0; k <= LAYERS+1; k++) {
                Temperature_last[i][j][k] = 0.0;
            }
        }
    }

    // Local boundry condition endpoints
    // Side Layers for each PE
    for(i = 0; i <= ROWS+1; i++){
        for(j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j][LAYERS+1] = 100;
        }
        for(k = 0; k <= LAYERS+1; k++) {
            Temperature_last[i][COLUMNS+1][k] = 100;
        }
    }


    // Bottom Layer
    if (myid == npes - 1){
        for(j = 0; j <= COLUMNS+1; j++){
            for(k = 0; k <= LAYERS+1; k++){
                Temperature_last[ROWS+1][j][k] = 100;
            }
        }
    }
}


// only called by last PE
void track_progress(int iteration, double dtt) {

    int i;

    printf("---------- Iteration number: %d , %f------------\n", iteration,dtt);

    // output global coordinates so user doesn't have to understand decomposition
/*     for(i = 5; i >= 0; i--) {
      printf("[%d,%d,%d]: %5.2f  ", ROWS_GLOBAL-i, COLUMNS-i, LAYERS - i, Temperature[ROWS-i][COLUMNS-i][LAYERS - i]);
    } */
    printf("\n");
}
