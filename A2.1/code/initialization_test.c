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
  
      int *loc_global_index;
      int *rank = (int*) malloc(sizeof(int)*(*nintcf + 1));
      int nintci_loc, nintcf_loc, nextci_loc, nextcf_loc;
      int r;
    
      for (r=0; r<4; r++)
      {
        
        allread_calc_global_idx(&loc_global_index, &nintci_loc, &nintcf_loc, &nextci_loc,
    			    &nextcf_loc, part_type, read_type, 4, r,
    			    *nintci, *nintcf, *nextci,
    			    *nextcf, *lcc, *elems, *points_count);
    
        for (i=nintci_loc; i <= nintcf_loc; i++) {
          rank[loc_global_index[nintci_loc + i]] = r;
        }
        
        free(loc_global_index); 
      }
  
  
  // For oneread
  
//   int **loc_global_index;
//   int *rank = (int*) malloc(sizeof(int)*(*nintcf + 1));
//   int *nintci_loc, *nintcf_loc, *nextci_loc, *nextcf_loc;
//   int r;
//   
//   oneread_calc_global_idx(&loc_global_index, &nintci_loc, &nintcf_loc, &nextci_loc,
// 			  &nextcf_loc, part_type, read_type, 5,
// 			  *nintci, *nintcf, *nextci,
// 			  *nextcf, *lcc, *elems, *points_count);
//   
//   for (r=0;r<5;r++)
//   {
//     printf("%d\t%d\n", r, nintcf_loc[r]);
//     for (i=nintci_loc[r]; i <= nintcf_loc[r]; i++) {
// //       rank[loc_global_index[r][nintci_loc[r] + i]] = r;
//     }
//   }
//   
//   //TODO: free it as a 2D array
//   //   for (r=0; r<4; r++)
//   //   {
//   //    free(loc_global_index[r]); 
//   //   }
//   //   free(loc_global_index); 
//   free(nintcf_loc);
//   free(nintci_loc);
//   free(nextci_loc);
//   free(nextcf_loc);
  
  
  
     vtk_write_unstr_grid_header("a", "b.vtk", *nintci, *nintcf, *points_count, *points, *elems);
     vtk_append_integer("b.vtk", "rank", *nintci, *nintcf, rank);
    
    free(rank);
  
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
		   
		   