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
void METIS_Partitioning(int **epart_ret, int *ne_ret, int nprocs, int *elems, int nintci, int nintcf, int points_count, int dual){
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  
  int NC, i;
  idx_t ncommon, nparts;
  idx_t *objval;
  idx_t ne, nn;
  
  ne = nintcf - nintci + 1;
  nn = points_count;
  ncommon = 4;
  nparts = nprocs;
  
  idx_t *epart;
  
  
  
  idx_t *eind;
  idx_t *eptr;
  
  //Unbelieveable thing happening - order matters!
  
   //eind = (idx_t *) calloc(sizeof(idx_t), (ne * 8));              
   //eptr = (idx_t *) calloc(sizeof(idx_t), (ne + 1));
  
  eptr = (idx_t *) calloc(sizeof(idx_t), (ne + 1));
  eind = (idx_t *) calloc(sizeof(idx_t), (8*ne));
  
  
  
  for (NC=0; NC<(ne+1); NC++) {
    eptr[NC] = 8*NC;
  }
  
  for (i=0; i<(8*ne); i++) {
    eind[i] = elems[i];
  }
  
  idx_t *npart;
  
  epart = (idx_t *) calloc(sizeof(idx_t), ne);
  
  npart = (idx_t *) calloc(sizeof(idx_t), nn);
  
  objval = (idx_t *) calloc(sizeof(idx_t), 1);
  
  if (dual == 1)
  {
    METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, NULL, objval, epart, npart);
  }
  else if (dual == 0)
  {
    METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL, &nparts, NULL, NULL, objval, epart, npart);
    
  }
  *epart_ret = (int *) calloc(sizeof(int), ne);
  *ne_ret = (nintcf - nintci + 1);
  
  for (i=0; i<ne; i++) {
    (*epart_ret)[i] = (int) (epart)[i];
  }
  
  
  free(epart);
  free(npart);
  free(eptr);
  free(eind);
  free(objval);
}

// For allread case - run on each process
void allread_calc_global_idx(int** local_global_index, int **global_local_index, int *nintci_loc, int *nintcf_loc, int *nextci_loc,
			     int *nextcf_loc, char *part_type, char*read_type, int nprocs, int myrank,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int *** send_lst, int **recv_cnt, int *** recv_lst )
{
  // Assume 0 value of nintci
  *nintci_loc = 0;
  int i, j, NC;
  int type, dual;
  determine_type(&type, &dual, part_type);  
  int el_count=0;
  int *epart;
  int ne;
 
  // If the classic mode chosen
  if (type == 0) {
    int start_int, stop_int, quotient_int, remainder_int, num_terms_int, num_terms_ext;
    ne= (nintcf - nintci + 1);

    epart = (int*) malloc( (ne)*sizeof(int) );
    
    // Calculation of the indices of the int cells for each process    
    num_terms_int = nintcf - nintci + 1;
    
    quotient_int  = num_terms_int  / nprocs;
    remainder_int = num_terms_int  % nprocs;
    
    start_int = (myrank + 0)*quotient_int  + MIN(myrank , remainder_int );
    stop_int = (myrank + 1)*quotient_int  + MIN(myrank+1, remainder_int );
    
    *nintcf_loc = stop_int - start_int - 1;
    *nextci_loc = stop_int - start_int;        
 
    // Calculation of the number of external cells belonging to a process
    num_terms_ext = 0;
    
    for (NC = start_int; NC < stop_int; NC++) {
      for (i=0; i<6; i++) {
	if (lcc[NC][i] >= nextci) {
	  (num_terms_ext)++;
	}
      }
    }
    
    *nextcf_loc = *nextci_loc + num_terms_ext - 1;    
    *local_global_index = (int*)malloc( (num_terms_int + num_terms_ext)*sizeof(int) );
    *global_local_index = (int*)malloc((nintcf+1)*sizeof(int));
       int k = 0;
       for (k=0; k<nintcf+1; k++){
           (*global_local_index)[k]=-1;
       }
    
    for (i=*nintci_loc; i <= *nintcf_loc; i++) {
      (*local_global_index)[i] = start_int + i;
      (*global_local_index)[start_int + i]=i;
    }
    
    // Calculation of the number of external cells belonging to a process
    j = *nextci_loc;
    
    for (NC = start_int; NC < stop_int; NC++) {
      for (i=0; i<6; i++) {
	if (lcc[NC][i] >= nextci) {
	  (*local_global_index)[j] = lcc[NC][i];
	  j++;
	}
      }
    }
    
    
    // Calculation of the epart for classic
    int rank;

    for (rank=0; rank<nprocs; rank++) {   
      start_int = (rank + 0)*quotient_int  + MIN(rank , remainder_int );
      stop_int = (rank + 1)*quotient_int  + MIN(rank+1, remainder_int );
      for (i=start_int; i<stop_int; i++) {
	epart[i] = rank;
      }
    }
    


        for (NC=0; NC<ne; NC++) {
            if (epart[NC] == myrank) {
    	    el_count++;
            }
        }
      
    
  } //if (type == 0)
  
  // If the metis mode chosen
  else if (type == 1) {

    METIS_Partitioning(&epart, &ne, nprocs, elems, nintci, nintcf, points_count, dual);    
    el_count=0;
    
    for (NC=0; NC<ne; NC++) {      
        if (epart[NC] == myrank) {
	    el_count++;
        }
    }
    
    *nintcf_loc = el_count-1;
    *nintci_loc = 0;
    int num_terms_ext;

    // Calculation of the number of external cells belonging to a process
    num_terms_ext = 0;
    
    for (NC = 0; NC < ne; NC++) {
      if (epart[NC] == myrank) {
	for (i=0; i<6; i++) {
	  if (lcc[NC][i] >= nextci) {
	    (num_terms_ext)++;
	  }
	}
      }
    }
    
   *local_global_index = (int*)malloc( (el_count + num_terms_ext)*sizeof(int) );
   /*allocate memory and initlaize global_local_index with -1*/
    *global_local_index = (int*)malloc((nintcf+1)*sizeof(int));
    int k = 0;
    for (k=0; k<nintcf+1; k++){
        (*global_local_index)[k]=-1;
    }
    
    *nextci_loc = el_count;
    *nextcf_loc = el_count + num_terms_ext - 1;
   
    i = 0;
    for (NC=0; NC<ne; NC++) {
      if (epart[NC] == myrank) {
	(*local_global_index)[i] =  NC;
	(*global_local_index)[NC] = i;
	i++;
      }
    }


    // Calculation of the number of external cells belonging to a process
    j = *nextci_loc;
    
    for (NC = 0; NC < ne; NC++) {
      if (epart[NC] == myrank) {
	for (i=0; i<6; i++) {
	  if (lcc[NC][i] >= nextci) {
	    (*local_global_index)[j] = lcc[NC][i];
	    j++;
	  }
	}
      }
    }
  }//else if (type == 1)

	/*neighbouring processor search */
    int *neighbour_proc_search;
    neighbour_proc_search = (int*)calloc(nprocs,sizeof(int));

    int global_index_temp;
    int current_neighbour;
    int proc_counter =0;

    for ( i =0; i< el_count; i++){
    		/*global index for the current local element */
    		global_index_temp = (*local_global_index)[i];
    		for (j =0; j<6; j++){
    			current_neighbour = lcc[global_index_temp][j];
    			if(current_neighbour > nintcf)
    				continue;
    			if(epart[current_neighbour]!=myrank){
    				neighbour_proc_search[epart[current_neighbour]] += 1;
    			}
    		}
    }

    for (i =0; i<nprocs; i++){
    	if(neighbour_proc_search[i]>0)
    		proc_counter +=1;
    }

    *nghb_cnt = proc_counter;
    *nghb_to_rank = (int*)malloc((*nghb_cnt)*sizeof(int));
    int *nghb_to_rank_reverse;
    nghb_to_rank_reverse = (int*)malloc(nprocs*sizeof(int));
    for (i =0; i<nprocs; i++){
    	nghb_to_rank_reverse[i]=-1;
    }
    int nghb_to_rank_counter = 0;
    for (i =0; i<nprocs; i++){
        	if(neighbour_proc_search[i]>0){
        		(*nghb_to_rank)[nghb_to_rank_counter] = i;
        		nghb_to_rank_reverse[i] = nghb_to_rank_counter;
        		nghb_to_rank_counter +=1 ;
        	}
        }

    /*quasi hash-table for the neighbour search*/
    /*very expensive!*/
	int **send_neighbour_search;
	send_neighbour_search = (int**)calloc((*nghb_cnt),sizeof(int*));
	for(i = 0; i < (*nghb_cnt); i++)
		send_neighbour_search[i] = (int*)calloc((nintcf+1),sizeof(int));

	int **recv_neighbour_search;
	recv_neighbour_search = (int**)calloc((*nghb_cnt),sizeof(int*));
	for(i = 0; i < (*nghb_cnt); i++)
		recv_neighbour_search[i] = (int*)calloc((nintcf+1),sizeof(int));

	/*allocate memory for send_cnt and revc_cnt*/
		*send_cnt = (int*)calloc((*nghb_cnt),sizeof(int));
		*recv_cnt = (int*)calloc((*nghb_cnt),sizeof(int));

    /*check total number of neighbours */

    for ( i =0; i< el_count; i++){
		/*global index for the current local element */        
		global_index_temp = (*local_global_index)[i];
		for (j =0; j<6; j++){
			current_neighbour = lcc[global_index_temp][j];

			if(current_neighbour > nintcf)
				continue;
			if(epart[current_neighbour]!= myrank){
				send_neighbour_search[nghb_to_rank_reverse[epart[current_neighbour]]][global_index_temp]+= 1;
				recv_neighbour_search[nghb_to_rank_reverse[epart[current_neighbour]]][current_neighbour]+= 1;
			}
		}
    }

	int destination =0;
	/*count the number of neighbors and the destination of sending*/
	for(proc_counter = 0; proc_counter < (*nghb_cnt); proc_counter++){
		for ( i = 0; i<nintcf+1; i++){
				if(send_neighbour_search[proc_counter][i] >0)
					(*send_cnt)[proc_counter]+=1;
				if(recv_neighbour_search[proc_counter][i] >0)
					(*recv_cnt)[proc_counter]+=1;
		}
	}

	/*memory allocation for nghb_to_rank, send_lst and recv_lst */
	if ( (*send_lst = (int**)malloc((*nghb_cnt) * sizeof(int*))) == NULL ) {
		fprintf(stderr, "malloc failed to allocate first dimension of send_lst");
	}
    for ( i = 0; i < (*nghb_cnt); i++ ) {
    	if ( ((*send_lst)[i] = (int *) malloc((*send_cnt)[i] * sizeof(int))) == NULL ) {
    		fprintf(stderr, "malloc failed to allocate second dimension of send_lst\n");
    	}
    }

	if ( (*recv_lst = (int**)malloc((*nghb_cnt)  * sizeof(int*))) == NULL ) {
		fprintf(stderr, "malloc failed to allocate first dimension of recv_lst");
	}
	for ( i = 0; i < (*nghb_cnt); i++ ) {
		if ( ((*recv_lst)[i] = (int *) malloc((*recv_cnt)[i] * sizeof(int))) == NULL ) {
			fprintf(stderr, "malloc failed to allocate second dimension of recv_lst\n");
		}
	}
	
	for(proc_counter = 0; proc_counter < (*nghb_cnt); proc_counter++){
		int send_counter = 0;
		int recv_counter = 0;
		for ( i = 0; i<nintcf+1; i++){
					if(send_neighbour_search[proc_counter][i] >0){
						(*send_lst)[proc_counter][send_counter]=i;
						send_counter ++;
					}
					if(recv_neighbour_search[proc_counter][i] >0){
						(*recv_lst)[proc_counter][recv_counter]=i;
						recv_counter++;
					}
		}
	}

	for (i =0; i<(*nghb_cnt); i++){
		free(send_neighbour_search[i]);
		free(recv_neighbour_search[i]);
	}

	free(send_neighbour_search);
	free(recv_neighbour_search);
    free(epart);
   //}//else if (type == 1)
}// allread_calc_global_idx



// For oneread case - run on 0 thread and distribute returned data structures in initialization.c
void oneread_calc_global_idx(int*** local_global_index, int ***global_local_index, int **nintci_loc, int **nintcf_loc, int **nextci_loc,
			     int **nextcf_loc, char *part_type, char*read_type, int nprocs,
			     int nintci, int nintcf, int nextci,int nextcf, int** lcc, int* elems, int points_count,
			     int **nghb_cnt, int ***nghb_to_rank, int *** send_cnt, int **** send_lst, int ***recv_cnt, int ****recv_lst)
{
  int i, j,NC;
  int type, dual;
  int rank;
  int destination;
  int send_counter = 0;
  int recv_counter = 0;
  int ne = (nintcf - nintci + 1);
  int *epart;
//  epart = (int*) malloc( (ne)*sizeof(int) );
  int* el_count;
  el_count = (int*)calloc( nprocs, sizeof(int) );
  //counter used in type = 1
  int *i_loc = (int*)calloc( nprocs, sizeof(int) );
  
  *nintci_loc = (int*)malloc( nprocs*sizeof(int) );
  *nintcf_loc = (int*)malloc( nprocs*sizeof(int) );
  *nextci_loc = (int*)malloc( nprocs*sizeof(int) );
  *nextcf_loc = (int*)malloc( nprocs*sizeof(int) );
  
  *local_global_index = (int**)malloc( nprocs*sizeof(int*) );
  *global_local_index = (int**)malloc( nprocs*sizeof(int*) );
  determine_type(&type, &dual, part_type); 
  // If the classic mode chosen
  if (type == 0) {
    int start_int, stop_int, quotient_int, remainder_int, num_terms_int, num_terms_ext;
    epart = (int*) malloc( (ne)*sizeof(int) );
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
      
      // Calculation of the epart for classic

      for (i=start_int; i<stop_int; i++) {
	epart[i] = rank;
      }
    
      
      // Calculation of the number of external cells belonging to a process
      num_terms_ext = 0;
      
      for (NC = start_int; NC < stop_int; NC++) {
	for (i=0; i<6; i++) {
	  if (lcc[NC][i] >= nextci) {
	    (num_terms_ext)++;
	  }
	}
      }
      
      (*nextcf_loc)[rank] = (*nextci_loc)[rank] + num_terms_ext - 1;     
      (*local_global_index)[rank] = (int*)malloc( (num_terms_int + num_terms_ext)*sizeof(int) );
      (*global_local_index)[rank] = (int*)malloc( (nintcf+1)*sizeof(int) );
      int k = 0;
          for (k=0; k<nintcf+1; k++){
              (*global_local_index)[rank][k]=-1;
          }
      for (i=(*nintci_loc)[rank]; i <= (*nintcf_loc)[rank]; i++) {
	(*local_global_index)[rank][i] = start_int + i;	
	(*global_local_index)[rank][start_int + i] = i;
      }      
      
      // Calculation of the number of external cells belonging to a process
      j = (*nextci_loc)[rank];
      
      for (NC = start_int; NC < stop_int; NC++) {
	for (i=0; i<6; i++) {
	  if (lcc[NC][i] >= nextci) {
	    (*local_global_index)[rank][j] = lcc[NC][i];
	    j++;
	  }
	}
      }
      
    }
    
    for (NC=0; NC<ne; NC++) {
          (el_count[epart[NC]])++;
        }
    
  }//if (type == 0)  

// If the metis mode chosen
  else if (type == 1) {
   // idx_t *epart;
    //idx_t ne;
    
    METIS_Partitioning(&epart, &ne, nprocs, elems, nintci, nintcf, points_count, dual);    


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
      num_terms_ext = 0;
      
      for (NC = 0; NC < ne; NC++) {
	if (epart[NC] == rank) {
	  for (i=0; i<6; i++) {
	    if (lcc[NC][i] >= nextci) {
	      (num_terms_ext)++;
	    }
	  }
	}
      }
      
      (*local_global_index)[rank] = (int*)malloc( (el_count[rank] + num_terms_ext)*sizeof(int) );
      (*global_local_index)[rank] = (int*)calloc( (nintcf+1),sizeof(int) );
     
      (*nextci_loc)[rank] = el_count[rank];
      (*nextcf_loc)[rank] = el_count[rank] + num_terms_ext - 1;
      
      // Calculation of the number of external cells belonging to a process
      j = (*nextci_loc)[rank];
      
      for (NC = 0; NC < ne; NC++) {
	if (epart[NC] == rank) {
	  for (i=0; i<6; i++) {
	    if (lcc[NC][i] >= nextci) {
	      (*local_global_index)[rank][j] = lcc[NC][i];
	      j++;
	    }
	  }
	}
      }
      
    }    

    
    for (NC=0; NC<ne; NC++) {
	(*local_global_index)[epart[NC]][i_loc[epart[NC]]] =  NC;
	(*global_local_index)[epart[NC]][NC] = i_loc[epart[NC]];
	(i_loc[epart[NC]])++; 
    }
    
  }// else if (type == 1)

    /*neighbouring processor search*/
    int **neighbour_proc_search;
        neighbour_proc_search = (int**)malloc(nprocs*sizeof(int*));
        for (i =0; i< nprocs; i++){
        	neighbour_proc_search[i] = (int*)calloc(nprocs,sizeof(int));
        }

        int global_index_temp;
        int current_neighbour;
        int proc_counter =0;
        int loop_counter = 0;
        int **nghb_to_rank_reverse;
        int nghb_to_rank_counter;

        /*memory allocation for the first dimension of array   int ***nghb_to_rank, int *** send_cnt, int **** send_lst, int ***recv_lst, int ****recv_lst*/
        *nghb_cnt = (int*)malloc(nprocs*sizeof(int));
        *nghb_to_rank = (int **)malloc(nprocs*sizeof(int*));
        nghb_to_rank_reverse = (int **)malloc(nprocs*sizeof(int*));
        *send_cnt = (int **)malloc(nprocs*sizeof(int*));
        *send_lst = (int ***)malloc(nprocs*sizeof(int**));
        *recv_cnt = (int **)malloc(nprocs*sizeof(int*));
        *recv_lst = (int ***)malloc(nprocs*sizeof(int**));

        for(loop_counter = 0; loop_counter < nprocs; loop_counter ++){
        	for ( i =0; i< el_count[loop_counter]; i++){
					/*global index for the current local element */
					global_index_temp = (*local_global_index)[loop_counter][i];
					for (j =0; j<6; j++){
						current_neighbour = lcc[global_index_temp][j];
						if(current_neighbour > nintcf)
							continue;
						if(epart[current_neighbour]!=loop_counter){
							neighbour_proc_search[loop_counter][epart[current_neighbour]] += 1;
						}
					}
			}
        }

        for(loop_counter = 0; loop_counter < nprocs; loop_counter ++){
        		 proc_counter = 0;
			for (i =0; i<nprocs; i++){
				if(neighbour_proc_search[loop_counter][i]>0)
					proc_counter +=1;
			}

			(*nghb_cnt)[loop_counter] = proc_counter;
			(*nghb_to_rank)[loop_counter] = (int*)malloc(((*nghb_cnt)[loop_counter])*sizeof(int));
			nghb_to_rank_reverse[loop_counter] = (int*)malloc(nprocs*sizeof(int));
			for (i =0; i<nprocs; i++){
				nghb_to_rank_reverse[loop_counter][i]=-1;
			}
			nghb_to_rank_counter = 0;
			for (i =0; i<nprocs; i++){
					if(neighbour_proc_search[loop_counter][i]>0){
						(*nghb_to_rank)[loop_counter][nghb_to_rank_counter] = i;
						nghb_to_rank_reverse[loop_counter][i] = nghb_to_rank_counter;
						nghb_to_rank_counter +=1 ;
					}
			}
        }

        for(loop_counter = 0; loop_counter < nprocs; loop_counter ++){

        /*quasi hash-table for the neighbour search*/
        /*very expensive!*/
    	int **send_neighbour_search;
    	send_neighbour_search = (int**)calloc((*nghb_cnt)[loop_counter],sizeof(int*));
    	for(i = 0; i < (*nghb_cnt)[loop_counter]; i++)
    		send_neighbour_search[i] = (int*)calloc((nintcf+1),sizeof(int));

    	int **recv_neighbour_search;
    	recv_neighbour_search = (int**)calloc((*nghb_cnt)[loop_counter],sizeof(int*));
    	for(i = 0; i < (*nghb_cnt)[loop_counter]; i++)
    		recv_neighbour_search[i] = (int*)calloc((nintcf+1),sizeof(int));

    	/*allocate memory for send_cnt and revc_cnt*/
    		(*send_cnt)[loop_counter] = (int*)calloc((*nghb_cnt)[loop_counter],sizeof(int));
    		(*recv_cnt)[loop_counter] = (int*)calloc((*nghb_cnt)[loop_counter],sizeof(int));

        /*check total number of neighbours */

        for ( i =0; i< el_count[loop_counter]; i++){
    		/*global index for the current local element */
    		global_index_temp = (*local_global_index)[loop_counter][i];
    		for (j =0; j<6; j++){
    			current_neighbour = lcc[global_index_temp][j];

    			if(current_neighbour > nintcf)
    				continue;
    			if(epart[current_neighbour]!= loop_counter){
    				send_neighbour_search[nghb_to_rank_reverse[loop_counter][epart[current_neighbour]]][global_index_temp]+= 1;
    				recv_neighbour_search[nghb_to_rank_reverse[loop_counter][epart[current_neighbour]]][current_neighbour]+= 1;
    			}
    		}
        }

    	destination =0;
    	/*count the number of neighbors and the destination of sending*/
    	for(proc_counter = 0; proc_counter < (*nghb_cnt)[loop_counter]; proc_counter++){
    		for ( i = 0; i<nintcf+1; i++){
    				if(send_neighbour_search[proc_counter][i] >0)
    					(*send_cnt)[loop_counter][proc_counter]+=1;
    				if(recv_neighbour_search[proc_counter][i] >0)
    					(*recv_cnt)[loop_counter][proc_counter]+=1;
    		}
    	}

    	/*memory allocation for nghb_to_rank, send_lst and recv_lst */
    	if ( ((*send_lst)[loop_counter] = (int**)malloc((*nghb_cnt)[loop_counter] * sizeof(int*))) == NULL ) {
    		fprintf(stderr, "malloc failed to allocate first dimension of send_lst");
    	}
        for ( i = 0; i < (*nghb_cnt)[loop_counter]; i++ ) {
        	if ( ((*send_lst)[loop_counter][i] = (int *) malloc((*send_cnt)[loop_counter][i] * sizeof(int))) == NULL ) {
        		fprintf(stderr, "malloc failed to allocate second dimension of send_lst\n");
        	}
        }

    	if ( ((*recv_lst)[loop_counter] = (int**)malloc((*nghb_cnt)[loop_counter] * sizeof(int*))) == NULL ) {
    		fprintf(stderr, "malloc failed to allocate first dimension of recv_lst");
    	}
    	for ( i = 0; i < (*nghb_cnt)[loop_counter]; i++ ) {
    		if ( ((*recv_lst)[loop_counter][i] = (int *) malloc((*recv_cnt)[loop_counter][i] * sizeof(int))) == NULL ) {
    			fprintf(stderr, "malloc failed to allocate second dimension of recv_lst\n");
    		}
    	}

    	for(proc_counter = 0; proc_counter < (*nghb_cnt)[loop_counter]; proc_counter++){
    		send_counter = 0;
    		recv_counter = 0;
    		for ( i = 0; i<nintcf+1; i++){
    					if(send_neighbour_search[proc_counter][i] >0){
    						(*send_lst)[loop_counter][proc_counter][send_counter]=i;
    						send_counter ++;
    					}
    					if(recv_neighbour_search[proc_counter][i] >0){
    						(*recv_lst)[loop_counter][proc_counter][recv_counter]=i;
    						recv_counter++;
    					}
    		}
    	}

    	for (i =0; i<(*nghb_cnt)[loop_counter]; i++){
    		free(send_neighbour_search[i]);
    		free(recv_neighbour_search[i]);
    	}
    	free(send_neighbour_search);
    	free(recv_neighbour_search);

        }//loop_counter




    free(i_loc);    
    free(el_count);    
    free(epart);        
  
}//oneread_calc_global_idx
