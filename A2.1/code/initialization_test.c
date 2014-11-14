/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdlib.h>
#include <stdio.h>
#include "util_read_files.h"
#include "initialization.h"
#include "domain_distribution.h"
#include "util_write_files.h"
#include "test_functions.h"

int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
		   int* nintci, int* nintcf, int* nextci,
		   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
		   double** bl, double** bh, double** bp, double** su, int* points_count,
		   int*** points, int** elems, double** var, double** cgup, double** oc,
		   double** cnorm, int** local_global_index) {
  /********** START INITIALIZATION **********/
  printf("INIT!\n");
  int i = 0;
  // read-in the input file
  int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
				 &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
				 &*points, &*elems);
  
  
  // ====================== For testing purposes ======================
  
  // For allread  
  
  //       int *loc_global_index;
  //       int *rank = (int*) malloc(sizeof(int)*(*nintcf + 1));
  //       int nintci_loc, nintcf_loc, nextci_loc, nextcf_loc;
  //       int r;
  //       
  //       int numproc = 3;
  //     
  //       for (r=0; r<numproc; r++)
  //       {
  //         
  //         allread_calc_global_idx(&loc_global_index, &nintci_loc, &nintcf_loc, &nextci_loc,
  //     			    &nextcf_loc, part_type, read_type, numproc, r,
  //     			    *nintci, *nintcf, *nextci,
  //     			    *nextcf, *lcc, *elems, *points_count);
  //     
  //         for (i=nintci_loc; i <= nintcf_loc; i++) {
  // 	  rank[loc_global_index[i - nintci_loc]] = r;
  //         }
  //         
  //         free(loc_global_index); 
  //       }
  //   
  
  // For oneread
  
  int **loc_global_index;
  int *rank = (int*) malloc(sizeof(int)*(*nintcf + 1));
  int *nintci_loc, *nintcf_loc, *nextci_loc, *nextcf_loc;
  int r;
  int numproc = 4;
  
  oneread_calc_global_idx(&loc_global_index, &nintci_loc, &nintcf_loc, &nextci_loc,
			  &nextcf_loc, part_type, read_type, numproc,
			  *nintci, *nintcf, *nextci,
			  *nextcf, *lcc, *elems, *points_count);
  
  for (r=0;r<numproc;r++)
  {
    printf("%d\t%d\n", r, nintcf_loc[r]);
    for (i=nintci_loc[r]; i <= nintcf_loc[r]; i++) {
      rank[(loc_global_index[r][i - nintci_loc[r]])] = r;
    }
  }
  
  
//     for (r=0; r<numproc; r++)
//   {
//     free(loc_global_index[r]); 
//   }
//   free(loc_global_index); 
//   free(nintcf_loc);
//   free(nintci_loc);
//   free(nextci_loc);
//   free(nextcf_loc);
  
  
  
  //      vtk_write_unstr_grid_header("a", "b.vtk", *nintci, *nintcf, *points_count, *points, *elems);
  //      vtk_append_integer("b.vtk", "rank", *nintci, *nintcf, rank);
  
  double *scalars = (double*) calloc(nextci_loc[2], sizeof(double));
  
  for (i=0; i<nextci_loc[2]; i++)
  {
    scalars[i] = 1.0;
  }
  
  test_distribution(file_in, "bb.vtk", loc_global_index[2], nextci_loc[2], scalars);
  free(scalars);
  free(rank);
  
  
  
    for (r=0; r<numproc; r++)
  {
    free(loc_global_index[r]); 
  }
  free(loc_global_index); 
  free(nintcf_loc);
  free(nintci_loc);
  free(nextci_loc);
  free(nextcf_loc);
  
  // ====================== For testing purposes ======================
  
  
  if ( f_status != 0 ) return f_status;
  
  *var = (double*) calloc(sizeof(double), (*nextcf + 1));
  *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));
  *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));
  
  // initialize the arrays
  for ( i = 0; i <= 10; i++ ) {
    (*cnorm)[i] = 1.0;
  }
  
  for ( i = (*nintci); i <= (*nintcf); i++ ) {
    (*var)[i] = 0.0;
  }
  
  for ( i = (*nintci); i <= (*nintcf); i++ ) {
    (*cgup)[i] = 1.0 / ((*bp)[i]);
  }
  
  for ( i = (*nextci); i <= (*nextcf); i++ ) {
    (*var)[i] = 0.0;
    (*cgup)[i] = 0.0;
    (*bs)[i] = 0.0;
    (*be)[i] = 0.0;
    (*bn)[i] = 0.0;
    (*bw)[i] = 0.0;
    (*bh)[i] = 0.0;
    (*bl)[i] = 0.0;
  }
  
  return 0;
		   }
		   
		   