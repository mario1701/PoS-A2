/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <mpi.h>
#include "util_write_files.h"

void finalization(char* file_in, int nprocs, int myrank, int total_iters, double residual_ratio,
                  int nintci, int nintcf, double* var, int* local_global_index, int* global_local_index) {

    if (nprocs == 1) {
      
    char file_out[100];
    sprintf(file_out, "%s_summary.out", file_in);
    //sprintf(file_out, "%s_summary.rank%d.out", file_in, myrank);
    
    int status = store_simulation_stats(file_in, file_out, nintci, nintcf, var, total_iters,
					residual_ratio);

    if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);
    
  }
  else {
  
    MPI_Request request;
  
    int global_nintci = 0;
    int global_num_internal_cells = 0;
    int global_nintcf = 0;
    int num_internal_cells = nintcf + 1;
    
    MPI_Reduce(&num_internal_cells, &global_num_internal_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    global_nintcf = global_num_internal_cells - 1;
    
    if (myrank > 0) {
      MPI_Isend(&nintcf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
      MPI_Isend(var, (nintcf+1), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
      MPI_Isend(local_global_index, (nintcf+1), MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
    }
  
     if (myrank == 0) {
       
       
       int i;
       double *global_var;
       int proc = 0;
       int nintcf_loc = 0;
       int ref_pos = nintcf+1;
       
       global_var = (double*) calloc(sizeof(double), (global_nintcf - global_nintci + 1));
       
       for (i = 0; i < nintcf + 1; i++) {
	 global_var[local_global_index[i]] = var[i];
       }
       
       for (proc = 1; proc < nprocs; proc++) {
	MPI_Recv(&nintcf_loc, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	double *var_loc = (double*) calloc(sizeof(double), (nintcf_loc + 1));
	int *local_global_index_loc = (int*) calloc(sizeof(int), (nintcf_loc + 1));
	
	//MPI_Recv(&(global_var[ref_pos]), (nintcf_loc+1), MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(var_loc, (nintcf_loc+1), MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(local_global_index_loc, (nintcf_loc+1), MPI_INT, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	for (i = 0; i < nintcf_loc + 1; i++) {
	  global_var[local_global_index_loc[i]] = var_loc[i];
	}
	
	free(var_loc);
	free(local_global_index_loc);
	
	//ref_pos += (nintcf_loc+1);
       }
    
      char file_out[100];
      sprintf(file_out, "%s_summary.out", file_in);
      //sprintf(file_out, "%s_summary.rank%d.out", file_in, myrank);
      

      int status = store_simulation_stats(file_in, file_out, global_nintci, global_nintcf, global_var, total_iters,
					  residual_ratio);
//       int status = store_simulation_stats(file_in, file_out, nintci, nintcf, var, total_iters,
// 					  residual_ratio);

      if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);
      
      free(global_var);
    
     }
  }

}

