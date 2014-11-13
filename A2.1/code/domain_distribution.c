/**
 * Domain distribution
 *
 * @date 10-Nov-2014
 * @author M. Bujny, Y. Dong
 */

#include <stdio.h>
#include "domain_distribution.h"
#include <malloc.h>
#include <metis.h>
#include <string.h>

void compute_boundary_start(int** boundary_direct_access, int *num_terms_ext, int nextcf, int nextci, int nintci_loc, int nintcf_loc, int **lcc)  {
  
  // Direct access table
  *boundary_direct_access = (int*)calloc( (nextcf - nextci + 1), sizeof(int) );
  int i, NC;
  
  for (NC = nintci_loc; NC <= nintcf_loc; NC++) {
    for (i=0; i<6; i++) {
      if (lcc[NC][i] >= nextci) {
	((*boundary_direct_access)[lcc[NC][i] - nextci])++;
      }
    }
  }
  
  *num_terms_ext = 0;
  
  for (NC=nextci; NC<=nextcf; NC++) {
    
    if ( (*boundary_direct_access)[NC - nextci] > 0) {
      (*num_terms_ext)++;
    }
    
  }
  
}

void compute_boundary_stop(int** boundary_direct_access, int *local_global_index, int nextcf, int nextci, int nextci_loc, int nextcf_loc, int **lcc) {
  
  int j=0;
  int i;
  
  for (i=nextci; i<=nextcf; i++) {
    if ( (*boundary_direct_access)[i - nextci] > 0) {
      local_global_index[nextci_loc + j] = i;
      //printf("%d\t%d\n",nextci_loc + j, local_global_index[nextci_loc + j]);
      j++;
    }
  }
  
  free(*boundary_direct_access);
  
}


void allread_calc_global_idx(int** local_global_index, int *nintci_loc, int *nintcf_loc, int *nextci_loc,
			     int *nextcf_loc, char *part_type, char*read_type, int nprocs, int myrank,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count) {
  
  printf("ALLREAD!\n");
  
  *nintci_loc = 0;
  int i, j, NC;
  int type, dual;
  
  if(strcmp(part_type, "classic")==0){
    type = 0;
    dual = 0;
  }
  
  if(strcmp(part_type, "dual")==0){
    type = 1;
    dual = 1;
  }
  
  if(strcmp(part_type, "nodal")==0){
    type = 1;
    dual = 0;
  }
  
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
    
    
    // Calculation of the number of external cells belonging to a process
    int *boundary_direct_access;
    compute_boundary_start(&boundary_direct_access, &num_terms_ext, nextcf, nextci, *nintci_loc, *nintcf_loc, lcc);
    
    
    *nextcf_loc = *nextci_loc + num_terms_ext - 1;
    
    *local_global_index = (int*)malloc( (num_terms_int + num_terms_ext)*sizeof(int) );
    
    
    for (i=*nintci_loc; i <= *nintcf_loc; i++) {
      (*local_global_index)[i] = start_int + i;
      
    }
    
    compute_boundary_stop(&boundary_direct_access, *local_global_index, nextcf, nextci, *nextci_loc, *nextcf_loc, lcc);
    
  } 
  
  else if (type == 1) {
    
    printf("METIS!\n");
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

printf("");
    
    if (dual == 1)
    {
      METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval, epart, npart);
 printf("ok!!!!!!!\n");  
    }
    else if (dual == 0)
    {
      METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL, &nparts, NULL, NULL, &objval, epart, npart);
    }
    
    int el_count=0;

    
    for (NC=0; NC<ne; NC++) {
      
      //  printf("\n%d\t%d", NC, epart[NC]);
      if (epart[NC] == myrank) {
	el_count++;
      }
      
    }
    
    *nintcf_loc = el_count-1;
    *nintci_loc = 0;
    
    int num_terms_ext;
    
    // Calculation of the number of external cells belonging to a process
    int *boundary_direct_access;
    compute_boundary_start(&boundary_direct_access, &num_terms_ext, nextcf, nextci, *nintci_loc, *nintcf_loc, lcc);
    
    printf("Elcount: %d\n", el_count);
    *local_global_index = (int*)malloc( (el_count + num_terms_ext)*sizeof(int) );
    

    *nextci_loc = el_count;
    *nextcf_loc = el_count + num_terms_ext - 1;
    
    
    i = 0;
    for (NC=0; NC<ne; NC++) {
      
      if (epart[NC] == myrank) {
	(*local_global_index)[i] =  NC;
	i++;
      }
      
    }
    
    compute_boundary_stop(&boundary_direct_access, *local_global_index, nextcf, nextci, *nextci_loc, *nextcf_loc, lcc);
    
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
			     
			     
