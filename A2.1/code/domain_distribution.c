/**
 * Domain distribution
 *
 * @date 10-Nov-2014
 * @author M. Bujny, Y. Dong
 */

#include <stdio.h>
#include "domain_distribution.h"
#include <malloc.h>
//#include "/usr/local/include/metis.h"
#include <metis.h>

void allread_calc_global_idx(int** local_global_index, int *nintci_loc, int *nintcf_loc, int *nextci_loc,
			     int *nextcf_loc, int type, int dual, int nprocs, int myrank,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count) {
  
  *nintci_loc = 0;
  int i, NC;
  
  if (type == 0) {
    int start_int, stop_int, quotient_int, remainder_int, num_terms_int;
    int start_ext, stop_ext, quotient_ext, remainder_ext, num_terms_ext;
    
    
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
    

    for (i=*nintci_loc; i <= *nintcf_loc; i++) {
      (*local_global_index)[i] = start_int + i;

    }
    
    //  printf("\nBoundary cells:\n");
    for (i=*nextci_loc; i <= *nextcf_loc; i++) {
      (*local_global_index)[i] = start_ext + i;

    }
    
  } 
  
  else if (type == 1) {
    
    
    idx_t options[METIS_NOPTIONS];
    //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
    
    METIS_SetDefaultOptions(options);
    
    idx_t ne, nn, ncommon, objval, nparts ;
    
    ne = (nintcf - nintci + 1);//(nintcf - nintci + 2);
    nn = points_count; //8*(nintcf - nintci + 1); 
    
    idx_t *eind = (idx_t*) malloc( 8*ne*sizeof(idx_t) );
    idx_t *eptr = (idx_t*) malloc( (ne+1)*sizeof(idx_t) );
    
    for (NC=0; NC<ne+1; NC++) {
      
      eptr[NC] = 8*NC;
      
    }
    
    for (i=0; i<8*ne; i++) {
      eind[i] = elems[i];
    }
    
    idx_t *epart;
    idx_t *npart;
    
    epart = (idx_t*) malloc( ne*sizeof(idx_t) );
    npart = (idx_t*) malloc( nn*sizeof(idx_t) );
    
    ncommon = 4;
    nparts = nprocs;
    
    METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval, epart, npart);
    
    int el_count=0;
    
    for (NC=0; NC<ne; NC++) {
      
      //  printf("\n%d\t%d", NC, epart[NC]);
      if (epart[NC] == myrank) {
	el_count++;
      }
      
    }
    
    printf("Elcount: %d\n", el_count);
    *local_global_index = (int*)malloc( (el_count)*sizeof(int) );
    
    
    
    *nintcf_loc = el_count-1;
    *nintci_loc = 0;
    
    i = 0;
    for (NC=0; NC<ne; NC++) {
      
      if (epart[NC] == myrank) {
	(*local_global_index)[i] =  NC;
	i++;
      }
      
    }
 
    free(epart);
    free(npart);
    free(eptr);
    free(eind);
    
    //free(local_global_index);
    
  }
  
  
  
  //free(local_global_index); 
			     }
			     
			     
			     /*
			      *			      n i*ntci_loc = malloc( nprocs*sizeof(int) );
			      *			      nintcf_loc = malloc( nprocs*sizeof(int) );
			      *			      nextci_loc = malloc( nprocs*sizeof(int) );
			      *			      nextcf_loc = malloc( nprocs*sizeof(int) );
			      */
			     