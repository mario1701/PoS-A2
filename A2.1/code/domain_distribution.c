/**
 * Domain distribution
 *
 * @date 10-Nov-2014
 * @author M. Bujny, Y. Dong
 */

#include <stdio.h>
#include "domain_distribution.h"
//#include <malloc.h>

void allread_calc_global_idx(int** local_global_index, int *nintci_loc, int *nintcf_loc, int *nextci_loc,
                   int *nextcf_loc, int type, int nprocs, int myrank,
                   int nintci, int nintcf, int nextci,
                   int nextcf, int** lcc) {
  
  *nintci_loc = 0;
  
  if (type == 0) {
    int start_int, stop_int, quotient_int, remainder_int, num_terms_int;
    int start_ext, stop_ext, quotient_ext, remainder_ext, num_terms_ext;
    int i;
  
    // Calculation of the indices of the int cells for each process
    
    num_terms_int = nintcf - nintci + 1;
  
    quotient_int  = num_terms_int  / nprocs;
    remainder_int = num_terms_int  % nprocs;

    start_int = (myrank + 0)*quotient_int  + MIN(myrank , remainder_int );
    stop_int = (myrank + 1)*quotient_int  + MIN(myrank+1, remainder_int );
    
    *nintcf_loc = stop_int - start_int - 1;
    *nextci_loc = stop_int - start_int;
    
    // Calculation of the indices of the int cells for each process
    
    num_terms_ext = nextcf - nextci + 1;
  
    quotient_ext  = num_terms_ext  / nprocs;
    remainder_ext = num_terms_ext  % nprocs;

    start_ext = *nextci_loc + (myrank + 0)*quotient_ext  + MIN(myrank , remainder_ext );
    stop_ext = *nextci_loc + (myrank + 1)*quotient_ext  + MIN(myrank+1, remainder_ext );
    
    *nextcf_loc = *nextci_loc + stop_ext - start_ext - 1;
    
    *local_global_index = (int*)malloc( ((stop_int - start_int) + (stop_ext - start_ext))*sizeof(int) );
    
    printf("\nInner cells:\n");
    for (i=*nintci_loc; i <= *nintcf_loc; i++) {
      (*local_global_index)[i] = start_int + i;
      printf("\n%d\t%d", i, (*local_global_index)[i]);
    }

    printf("\nBoundary cells:\n");
    for (i=*nextci_loc; i <= *nextcf_loc; i++) {
      (*local_global_index)[i] = start_ext + i;
      printf("\n%d\t%d", i, (*local_global_index)[i]);
    }
  
  } 
  //free(local_global_index); 
}

 
  /*
  nintci_loc = malloc( nprocs*sizeof(int) );
  nintcf_loc = malloc( nprocs*sizeof(int) );
  nextci_loc = malloc( nprocs*sizeof(int) );
  nextcf_loc = malloc( nprocs*sizeof(int) );
  */