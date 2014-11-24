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

// Start counting the external cells
// If an internal cell belonging to the particular process refers to an external cell, increase its counter
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
}//compute_boundary_start

// Append the external cells that were reffered at least once by internals cell to the mapping array
void compute_boundary_stop(int** boundary_direct_access, int *local_global_index, int nextcf, int nextci, int nextci_loc, int nextcf_loc, int **lcc) {
  int j=0;
  int i;
 
  for (i=nextci; i<=nextcf; i++) {
    if ( (*boundary_direct_access)[i - nextci] > 0) {
      local_global_index[nextci_loc + j] = i;
      j++;
    }
  }

  free(*boundary_direct_access);
}//compute_boundary_stop

// Deterine the type of data distribution
void determine_type(int *type, int *dual, char *part_type) {
  if(strcmp(part_type, "classic")==0){
    *type = 0;
    *dual = 0;
  }
  if(strcmp(part_type, "dual")==0){
    *type = 1;
    *dual = 1;
  }
  if(strcmp(part_type, "nodal")==0){
    *type = 1;
    *dual = 0;
  }
}

// Initialize METIS and devide the domain
void METIS_Partitioning(idx_t **epart, idx_t *ne, int nprocs, int *elems, int nintci, int nintcf, int points_count, int dual){
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    
    int NC, i;
    idx_t nn, ncommon, objval, nparts ;
    *ne = (nintcf - nintci + 1);
    nn = points_count; 
    
    idx_t *eind = (idx_t*) malloc( 8*(*ne)*sizeof(idx_t) );
    idx_t *eptr = (idx_t*) malloc( ((*ne)+1)*sizeof(idx_t) );
    
    for (NC=0; NC<(*ne)+1; NC++) {
      eptr[NC] = 8*NC;
    }
    
    for (i=0; i<8*(*ne); i++) {
      eind[i] = elems[i];
    }
    
    idx_t *npart;
    
    *epart = (idx_t*) malloc( (*ne)*sizeof(idx_t) );
    npart = (idx_t*) malloc( nn*sizeof(idx_t) );
    
    ncommon = 4;
    nparts = nprocs;

    if (dual == 1)
    {
      METIS_PartMeshDual(&(*ne), &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval, *epart, npart);
    }
    else if (dual == 0)
    {
      METIS_PartMeshNodal(&(*ne), &nn, eptr, eind, NULL, NULL, &nparts, NULL, NULL, &objval, *epart, npart);
    }
    
    free(npart);
    free(eptr);
    free(eind);
}

// For allread case - run on each process
void allread_calc_global_idx(int** local_global_index, int *nintci_loc, int *nintcf_loc, int *nextci_loc,
			     int *nextcf_loc, char *part_type, char*read_type, int nprocs, int myrank,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count) 
{
  // Assume 0 value of nintci
  *nintci_loc = 0;
  int i, j, NC;
  int type, dual;  
  determine_type(&type, &dual, part_type);  
 
  // If the classic mode chosen
  if (type == 0) {
    int start_int, stop_int, quotient_int, remainder_int, num_terms_int, num_terms_ext;    
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
  } //if (type == 0)
  
  // If the metis mode chosen
  else if (type == 1) {    
    idx_t *epart;
    idx_t ne;
    
    METIS_Partitioning(&epart, &ne, nprocs, elems, nintci, nintcf, points_count, dual);    
    int el_count=0;
    
    for (NC=0; NC<ne; NC++) {      
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
   }//else if (type == 1)
}// allread_calc_global_idx


// For oneread case - run on 0 thread and distribute returned data structures in initialization.c
void oneread_calc_global_idx(int*** local_global_index, int **nintci_loc, int **nintcf_loc, int **nextci_loc,
			     int **nextcf_loc, char *part_type, char*read_type, int nprocs,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count) 
{
  int i, NC;
  int type, dual;
  int rank;
  
  *nintci_loc = (int*)malloc( nprocs*sizeof(int) );
  *nintcf_loc = (int*)malloc( nprocs*sizeof(int) );
  *nextci_loc = (int*)malloc( nprocs*sizeof(int) );
  *nextcf_loc = (int*)malloc( nprocs*sizeof(int) );
  
  *local_global_index = (int**)malloc( nprocs*sizeof(int*) );  
  determine_type(&type, &dual, part_type); 
  // If the classic mode chosen
  if (type == 0) {
    int start_int, stop_int, quotient_int, remainder_int, num_terms_int, num_terms_ext;
    // Calculation of the indices of the int cells for each process
    num_terms_int = nintcf - nintci + 1;
    quotient_int  = num_terms_int  / nprocs;
    remainder_int = num_terms_int  % nprocs;

    for (rank=0; rank<nprocs; rank++) {
      // Assume 0 value of nintci
      (*nintci_loc)[rank] = 0;     
      start_int = (rank + 0)*quotient_int  + MIN(rank , remainder_int );
      stop_int = (rank + 1)*quotient_int  + MIN(rank+1, remainder_int );     
      (*nintcf_loc)[rank] = stop_int - start_int - 1;
      (*nextci_loc)[rank] = stop_int - start_int;
      // Calculation of the number of external cells belonging to a process
      int *boundary_direct_access;
      compute_boundary_start(&boundary_direct_access, &num_terms_ext, nextcf, nextci, (*nintci_loc)[rank], (*nintcf_loc)[rank], lcc);
      (*nextcf_loc)[rank] = (*nextci_loc)[rank] + num_terms_ext - 1;     
      (*local_global_index)[rank] = (int*)malloc( (num_terms_int + num_terms_ext)*sizeof(int) );

      for (i=(*nintci_loc)[rank]; i <= (*nintcf_loc)[rank]; i++) {
	(*local_global_index)[rank][i] = start_int + i;	
      }      
      compute_boundary_stop(&boundary_direct_access, (*local_global_index)[rank], nextcf, nextci, (*nextci_loc)[rank], (*nextcf_loc)[rank], lcc);
    }    
  }//if (type == 0)  

// If the metis mode chosen
  else if (type == 1) {
    idx_t *epart;
    idx_t ne;
    
    METIS_Partitioning(&epart, &ne, nprocs, elems, nintci, nintcf, points_count, dual);    
    int* el_count;    
    el_count = (int*)calloc( nprocs, sizeof(int) );

    for (NC=0; NC<ne; NC++) {      
      (el_count[epart[NC]])++;      
    }
    
    for (rank=0; rank<nprocs; rank++) {     
      (*nintcf_loc)[rank] = el_count[rank]-1;
      (*nintci_loc)[rank] = 0;
    }    
    
    for (rank=0; rank<nprocs; rank++) {
      int num_terms_ext;      
      // Calculation of the number of external cells belonging to a process
      int *boundary_direct_access;
      compute_boundary_start(&boundary_direct_access, &num_terms_ext, nextcf, nextci, (*nintci_loc)[rank], (*nintcf_loc)[rank], lcc);      
      (*local_global_index)[rank] = (int*)malloc( (el_count[rank] + num_terms_ext)*sizeof(int) );
     
      (*nextci_loc)[rank] = el_count[rank];
      (*nextcf_loc)[rank] = el_count[rank] + num_terms_ext - 1;
      
      compute_boundary_stop(&boundary_direct_access, (*local_global_index)[rank], nextcf, nextci, (*nextci_loc)[rank], (*nextcf_loc)[rank], lcc);      
    }    
        int* i_loc = (int*)calloc( nprocs, sizeof(int) );
    
    for (NC=0; NC<ne; NC++) {
	(*local_global_index)[epart[NC]][i_loc[epart[NC]]] =  NC;
	(i_loc[epart[NC]])++; 
    }
    
    free(i_loc);    
    free(el_count);    
    free(epart);        
  }//else if (type == 1) 
  
}//oneread_calc_global_idx
