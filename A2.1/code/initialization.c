/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

//#define PAPI

#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include "util_read_files.h"
#include <string.h>
#include "initialization.h"
#include <malloc.h>
#include "test_functions.h"
#include "domain_distribution.h"

#ifdef PAPI
#include <papi.h>
#endif
void write_vtk(char *file_in, char *scalars_name, int *local_global_index, int num_internal_cells, double *scalars, char *part_type, int myrank) ;

int memoryallocation(int ***LCC_local, double **bs_local, double **be_local, double **bn_local, double **bw_local, double **bh_local, double **bl_local, double **bp_local, double **su_local,int num_internal_cells, int num_cells, int *nintcf, int points_count, double **var, double **cgup, double **oc, double **cnorm);


int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
		   int* nintci, int* nintcf, int* nextci,
		   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
		   double** bl, double** bh, double** bp, double** su, int* points_count,
		   int*** points, int** elems, double** var, double** cgup, double** oc,
		   double** cnorm, int** local_global_index) 
{
  /********** START INITIALIZATION **********/
  
  #ifdef PAPI
  /*PAPI initialization*/
  float rtime, ptime, mflops;
  long long flpops;
 	void handle_error (int retval)
		{
		     printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
		  exit(1);
		}  
    if(PAPI_flops( &rtime, &ptime, &flpops,  &mflops ) != PAPI_OK) handle_error(1);
    /*decide key type*/
int input_key, part_key, read_key;

if(strcmp(read_type, "oneread") == 0){
read_key = 1;
}
if (strcmp(read_type, "allread") == 0){
read_key = 2;
}
if(strcmp(part_type, "classic") == 0){
part_key = 1;
}
if(strcmp(part_type, "dual") == 0){
part_key = 2;
}
if(strcmp(part_type, "nodal") == 0){
part_key = 3;
}
if((strcmp(file_in, "tjunc.geo.bin") ==0) ||(strcmp(file_in, "tjunc.geo.dat") ==0)  ){
input_key = 1;
}
if((strcmp(file_in, "drall.geo.bin") ==0) ||(strcmp(file_in, "drall.geo.dat") ==0)  ){
input_key = 2;
}
if((strcmp(file_in, "pent.geo.bin") ==0) ||(strcmp(file_in, "pent.geo.dat") ==0)  ){
input_key = 3;
}
if((strcmp(file_in, "cojack.geo.bin") ==0) ||(strcmp(file_in, "cojack.geo.dat") ==0)  ){
input_key = 4;
}
    #endif
  
  int i = 0;
  int j = 0;
 
  /***************read-in the input file*********************/
  if (strcmp(read_type, "oneread") == 0 ) {
    if(0 == myrank){
      int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
				     &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
				     &*points, &*elems);
      if ( f_status != 0 ) return f_status;
    }                               
  }
 
  if (strcmp(read_type, "allread") == 0){
    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
				   &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
				   &*points, &*elems);
    if ( f_status != 0 ) return f_status;
  }
 
  /************************data distribution*************************/
  
  
  /*local property*/
  int num_cells;
  int num_internal_cells;
  int Nintci_loc, Nintcf_loc, Nextci_loc, Nextcf_loc;/*local index for the beginning and ending of internal/external cells*/
  
  /*initialize the array which will only be contained locally*/ 
  int **LCC_local; 
  double *bs_local, *be_local, *bn_local, *bw_local, *bh_local, *bl_local;
  double *bp_local; 
  double *su_local; 
  //int local_global_points_count= *points_count;
  
  if(strcmp(read_type, "oneread") == 0){
     MPI_Status Status[6];
    int ** local_global_index_array;
    int *nintci_loc_array, *nintcf_loc_array, *nextci_loc_array, *nextcf_loc_array, *length_loc_index_array;
    int local_global_nintcf;
    int length_loc_index;
    
    if (0 == myrank){
      oneread_calc_global_idx(&local_global_index_array, &nintci_loc_array, &nintcf_loc_array, &nextci_loc_array,
			      &nextcf_loc_array, part_type, read_type, nprocs,
			      *nintci, *nintcf, *nextci,
			      *nextcf, *lcc, *elems, *points_count);
      int dest;
      length_loc_index_array = (int*) malloc(nprocs*sizeof(int));
      
      for (i=0; i<nprocs; i++) {
	length_loc_index_array[i] = nextcf_loc_array[i] - nintci_loc_array[i] + 1;
      }
      /*special treatment for processor 0 since there is no MPI_Send for that*/
      length_loc_index = length_loc_index_array[0];
      
      Nintci_loc = nintci_loc_array[0];
      Nintcf_loc = nintcf_loc_array[0];
      Nextci_loc = nextci_loc_array[0];
      Nextcf_loc = nextcf_loc_array[0];
      
      (*local_global_index) = (int*)malloc(length_loc_index*sizeof(int)); 
   
      for (i=0; i<length_loc_index; i++) {
	(*local_global_index)[i] = local_global_index_array[0][i];
      }
    
      for (dest = 1; dest<nprocs; dest++){
	MPI_Send(&(length_loc_index_array[dest]), 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
	MPI_Send(&(nintci_loc_array[dest]), 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
	MPI_Send(&(nintcf_loc_array[dest]), 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
	MPI_Send(&(nextci_loc_array[dest]), 1, MPI_INT, dest, 3, MPI_COMM_WORLD);
	MPI_Send(&(nextcf_loc_array[dest]), 1, MPI_INT, dest, 4, MPI_COMM_WORLD);
	MPI_Send(local_global_index_array[dest], length_loc_index_array[dest], MPI_INT, dest, 5, MPI_COMM_WORLD);
	MPI_Send(nintcf, 1, MPI_INT, dest, 6, MPI_COMM_WORLD);
	MPI_Send(points_count, 1, MPI_INT, dest, 7, MPI_COMM_WORLD);
      }
    }//if (0 == myrank)
    
    if (myrank>0){
      MPI_Recv(&length_loc_index,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      if ((*local_global_index = (int*)malloc(length_loc_index*sizeof(int)))== NULL){
	fprintf(stderr, "malloc(local_global_index) failed\n");
	return -1;
      }
      MPI_Recv(&Nintci_loc,1,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&Nintcf_loc,1,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&Nextci_loc,1,MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&Nextcf_loc,1,MPI_INT,0,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv((*local_global_index),length_loc_index,MPI_INT,0,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&local_global_nintcf, 1, MPI_INT, 0, 6, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&(*points_count), 1, MPI_INT, 0, 7, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }//if (myrank>0)
    
    /*Data Transfer*/
    if (myrank==0){
      int dest;
      for (dest = 1; dest<nprocs; dest++){
	num_cells = nextcf_loc_array[dest] - nintci_loc_array[dest] +1;
	num_internal_cells = nintcf_loc_array[dest] - nintci_loc_array[dest] +1;
	
	memoryallocation(&LCC_local, &bs_local, &be_local, &bn_local, &bw_local, &bh_local, &bl_local, &bp_local, &su_local, num_internal_cells, num_cells, &*nintcf, *points_count, &*var, &*cgup, &*oc, &*cnorm);
	
	for (i =nintci_loc_array[dest]; i <=nextcf_loc_array[dest]; i++){
	  bs_local[i] = (*bs)[local_global_index_array[dest][i]];
	  be_local[i] = (*be)[local_global_index_array[dest][i]];
	  bn_local[i] = (*bn)[local_global_index_array[dest][i]];
	  bw_local[i] = (*bw)[local_global_index_array[dest][i]];
	  bh_local[i] = (*bh)[local_global_index_array[dest][i]];
	  bl_local[i] = (*bl)[local_global_index_array[dest][i]];
	  bp_local[i] = (*bp)[local_global_index_array[dest][i]];
	  su_local[i] = (*su)[local_global_index_array[dest][i]];
	}

	MPI_Send(bs_local, num_cells, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
	MPI_Send(be_local, num_cells, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
	MPI_Send(bn_local, num_cells, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
	MPI_Send(bw_local, num_cells, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
	MPI_Send(bh_local, num_cells, MPI_DOUBLE, dest, 4, MPI_COMM_WORLD);
	MPI_Send(bl_local, num_cells, MPI_DOUBLE, dest, 5, MPI_COMM_WORLD);
	MPI_Send(bp_local, num_cells, MPI_DOUBLE, dest, 6, MPI_COMM_WORLD);
	MPI_Send(su_local, num_cells, MPI_DOUBLE, dest, 7, MPI_COMM_WORLD);
	
	for (i =nintci_loc_array[dest]; i <=nintcf_loc_array[dest]; i++){
	  for (j = 0; j<6; j++){
	    LCC_local[i][j] = (*lcc)[local_global_index_array[dest][i]][j];
	  }
	}

	for (i =nintci_loc_array[dest]; i <=nintcf_loc_array[dest]; i++){
	  MPI_Send(LCC_local[i], 6, MPI_INT, dest, i, MPI_COMM_WORLD);
	}

	for (i=0; i<num_internal_cells; i++) {
	  free(LCC_local[i]);
	}
	
	free(LCC_local);
	free(bs_local);
	free(be_local);
	free(bn_local); 
	free(bw_local);
	free(bh_local);
	free(bl_local);
	free(bp_local);
	free(su_local); 
	free(*var);
	free(*cgup);
	free(*oc);
	free(*cnorm);
      } //for (dest = 1; dest<nprocs; dest++)   
    }//if (myrank==0)
 
    if (myrank>0){
      num_cells = Nextcf_loc - Nintci_loc +1;
      num_internal_cells = Nintcf_loc  - Nintci_loc +1;
      
      memoryallocation(lcc, bs, be, bn, bw, bh, bl, bp, su, num_internal_cells, num_cells, &local_global_nintcf, (*points_count), &*var, &*cgup, &*oc, &*cnorm);
     
      MPI_Recv(*bs, num_cells, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(*be, num_cells, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(*bn, num_cells, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(*bw, num_cells, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(*bh, num_cells, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(*bl, num_cells, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(*bp, num_cells, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(*su, num_cells, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i =Nintci_loc; i <=Nintcf_loc; i++){
	MPI_Recv((*lcc)[i], 6, MPI_INT, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }//if (myrank>0)
   
    if (myrank==0){
      int dest = 0;
      num_cells = Nextcf_loc - Nintci_loc +1;
      num_internal_cells = Nintcf_loc  - Nintci_loc +1;
 
      memoryallocation(&LCC_local, &bs_local, &be_local, &bn_local, &bw_local, &bh_local, &bl_local, &bp_local, &su_local, num_internal_cells, num_cells, &*nintcf, *points_count, &*var, &*cgup, &*oc, &*cnorm);

      for (i =nintci_loc_array[dest]; i <=nextcf_loc_array[dest]; i++){
	bs_local[i] = (*bs)[local_global_index_array[dest][i]];
	be_local[i] = (*be)[local_global_index_array[dest][i]];
	bn_local[i] = (*bn)[local_global_index_array[dest][i]];
	bw_local[i] = (*bw)[local_global_index_array[dest][i]];
	bh_local[i] = (*bh)[local_global_index_array[dest][i]];
	bl_local[i] = (*bl)[local_global_index_array[dest][i]];
	bp_local[i] = (*bp)[local_global_index_array[dest][i]];
	su_local[i] = (*su)[local_global_index_array[dest][i]];
      }
      for (i =nintci_loc_array[dest]; i <=nintcf_loc_array[dest]; i++){
	for (j = 0; j<6; j++){
	  LCC_local[i][j] = (*lcc)[local_global_index_array[dest][i]][j];
        }
      }
 
      for (i=0; i<num_internal_cells; i++) {
	free((*lcc)[i]);
      }
      
      free(*lcc);
      free(*bs);
      free(*be);
      free(*bn); 
      free(*bw);
      free(*bh);
      free(*bl);
      free(*bp);
      free(*su); 
    
      *lcc = LCC_local;
      *bs = bs_local;
      *be = be_local;
      *bn = bn_local;
      *bw = bw_local;
      *bh = bh_local;
      *bl = bl_local;
      *bp = bp_local;
      *su = su_local;

      for (i=0; i<nprocs; i++) {
	free(local_global_index_array[i]);
      }
      free(local_global_index_array);
      free(length_loc_index_array);
      free(nintci_loc_array);
      free(nintcf_loc_array);
      free(nextci_loc_array);
      free(nextcf_loc_array); 
    }//if (myrank==0)

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
      (*cnorm)[i] = 1.0;
    }
    for ( i = Nintci_loc; i <= Nintcf_loc; i++ ) {
      (*var)[i] = 0.0;
    }
    for ( i = Nintci_loc; i <= Nintcf_loc; i++ ) {
      (*cgup)[i] = 1.0 / ((*bp)[i]);
    }
    for ( i = Nextci_loc; i <= Nextcf_loc; i++ ) {
      (*var)[i] = 0.0;
      (*cgup)[i] = 0.0;
      (*bs)[i] = 0.0;
      (*be)[i] = 0.0;
      (*bn)[i] = 0.0;
      (*bw)[i] = 0.0;
      (*bh)[i] = 0.0;
      (*bl)[i] = 0.0;
    }

#ifdef PAPI
	if(PAPI_flops( &rtime, &ptime, &flpops, &mflops ) != PAPI_OK) handle_error(1);
	write_pstats_exectime(input_key, part_key, read_key, myrank, ptime);
	write_pstats_partition(input_key, part_key, myrank, num_internal_cells, (Nextcf_loc - Nextci_loc +1) );
#endif

    write_vtk(file_in, "CGUP", *local_global_index, num_internal_cells, *cgup, part_type, myrank);
    write_vtk(file_in, "SU", *local_global_index, num_internal_cells, *su, part_type, myrank);

/*allocate memory for elems and points for processor >0 */
if ( (*elems = (int*) malloc((Nintcf_loc + 1) * 8 * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate elems");
        return -1;
    }
if ( (*points = (int **) calloc((*points_count), sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }

    for ( i = 0; i < (*points_count); i++ ) {
        if ( ((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    }


  }//if(strcmp(read_type, "oneread") == 0)
  
  
  
  if(strcmp(read_type, "allread") == 0){
    
    allread_calc_global_idx( &(*local_global_index), &Nintci_loc, &Nintcf_loc, &Nextci_loc, &Nextcf_loc, part_type, read_type, nprocs, myrank,*nintci, *nintcf, *nextci, *nextcf, *lcc, *elems, *points_count);     
    
    num_cells = Nextcf_loc - Nintci_loc +1;
    num_internal_cells = Nintcf_loc - Nintci_loc +1;
    
    memoryallocation(&LCC_local, &bs_local, &be_local, &bn_local, &bw_local, &bh_local, &bl_local, &bp_local, &su_local, num_internal_cells, num_cells, &*nintcf, *points_count, &*var, &*cgup, &*oc, &*cnorm);

    /*read LCC for LCC_local*/
    for (i =Nintci_loc; i <=Nintcf_loc; i++){
      for (j = 0; j<6; j++){
	LCC_local[i][j] = (*lcc)[(*local_global_index)[i]][j];
      }
    }
    /*read arrays*/
    for (i =Nintci_loc; i <=Nextcf_loc; i++){
      bs_local[i] = (*bs)[(*local_global_index)[i]];
      be_local[i] = (*be)[(*local_global_index)[i]];
      bn_local[i] = (*bn)[(*local_global_index)[i]];
      bw_local[i] = (*bw)[(*local_global_index)[i]];
      bh_local[i] = (*bh)[(*local_global_index)[i]];
      bl_local[i] = (*bl)[(*local_global_index)[i]];
      bp_local[i] = (*bp)[(*local_global_index)[i]];
      su_local[i] = (*su)[(*local_global_index)[i]];
    }
    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
      (*cnorm)[i] = 1.0;
    }
    for ( i = Nintci_loc; i <= Nintcf_loc; i++ ) {
      (*var)[i] = 0.0;
    }
    for ( i = Nintci_loc; i <= Nintcf_loc; i++ ) {
      (*cgup)[i] = 1.0 / (bp_local[i]);
    }
    for ( i = Nextci_loc; i <= Nextcf_loc; i++ ) {
      (*var)[i] = 0.0;
      (*cgup)[i] = 0.0;
      (*bs)[i] = 0.0;
      (*be)[i] = 0.0;
      (*bn)[i] = 0.0;
      (*bw)[i] = 0.0;
      (*bh)[i] = 0.0;
      (*bl)[i] = 0.0;
    }     
    
    /*exchange the memory name for local and global and free the global one*/
    int **LCC_local_temp=NULL; 
    double *bs_local_temp=NULL;
    double *be_local_temp = NULL;
    double *bn_local_temp = NULL;
    double *bw_local_temp = NULL; 
    double *bh_local_temp = NULL;
    double *bl_local_temp = NULL;
    double *bp_local_temp = NULL; 
    double *su_local_temp = NULL;
   
    LCC_local_temp = LCC_local; LCC_local = *lcc;  *lcc = LCC_local_temp;
    bs_local_temp = bs_local; bs_local = *bs;  *bs = bs_local_temp;
    be_local_temp = be_local; be_local = *be;  *be = be_local_temp;
    bn_local_temp = bn_local; bn_local = *bn;  *bn = bn_local_temp;
    bw_local_temp = bw_local; bw_local = *bw;  *bw = bw_local_temp;
    bh_local_temp = bh_local; bh_local = *bh;  *bh = bh_local_temp;
    bl_local_temp = bl_local; bl_local = *bl;  *bl = bl_local_temp;
    bp_local_temp = bp_local; bp_local = *bp;  *bp = bp_local_temp;
    su_local_temp = su_local; su_local = *su;  *su = su_local_temp;
 
#ifdef PAPI
	if(PAPI_flops( &rtime, &ptime, &flpops, &mflops ) != PAPI_OK) handle_error(1);
	write_pstats_exectime(input_key, part_key, read_key, myrank, ptime);
	write_pstats_partition(input_key, part_key, myrank, num_internal_cells, (Nextcf_loc - Nextci_loc +1) );
#endif
    
    double *ranks = (double*) calloc(num_internal_cells, sizeof(double));
    for (i=0; i<num_internal_cells; i++)
    {
      ranks[i] = 1.0;
    }
    write_vtk(file_in, "ranks", *local_global_index, num_internal_cells, ranks, part_type, myrank);
    write_vtk(file_in, "CGUP", *local_global_index, num_internal_cells, *cgup, part_type, myrank);
    write_vtk(file_in, "SU", *local_global_index, num_internal_cells, *su, part_type, myrank);
    
    free(ranks);

    free(su_local);
    free(bp_local);
    free(bh_local);
    free(bl_local);
    free(bw_local);
    free(bn_local);
    free(be_local);
    free(bs_local);
    for ( i = 0; i < (*nintcf + 1); i++ ) {
      free(LCC_local[i]);
    }
    free(LCC_local);  
  }//if(strcmp(read_type, "allread") == 0)
  
  //*points_count = local_global_points_count;
  
  *nintci = Nintci_loc;
  *nintcf = Nintcf_loc;
  *nextci = Nextci_loc;
  *nextcf = Nextcf_loc;

  return 0;
}

int memoryallocation(int ***LCC_local, double **bs_local, double **be_local, double **bn_local, double **bw_local, double **bh_local, double **bl_local, double **bp_local, double **su_local, /*int points_count_local,int ***points_local, int **elems_local,*/ int num_internal_cells, int num_cells, int *nintcf, int points_count, double **var, double **cgup, double **oc, double **cnorm)
{
  // allocating LCC_local
  int i =0;  
  if ( (*LCC_local = (int**) malloc((num_internal_cells) * sizeof(int*))) == NULL ) {
    fprintf(stderr, "malloc failed to allocate first dimension of LCC");
    return -1;
  }
  
  for ( i = 0; i < num_internal_cells; i++ ) {
    if ( ((*LCC_local)[i] = (int*) malloc(6 * sizeof(int))) == NULL ) {
      fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
      return -1;
    }
  }
  
  // allocate other arrays
  if ( (*bs_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(BS) failed\n");
    return -1;
  }
  
  if ( (*be_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(BE) failed\n");
    return -1;
  }
  
  if ( (*bn_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(BN) failed\n");
    return -1;
  }
  
  if ( (*bw_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(BW) failed\n");
    return -1;
  }
  
  if ( (*bh_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(BL) failed\n");
    return -1;
  }
  
  if ( (*bl_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(BH) failed\n");
    return -1;
  }
  
  if ( (*bp_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(BP) failed\n");
    return -1;
  }
  
  if ( (*su_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
    fprintf(stderr, "malloc(SU) failed\n");
    return -1;
  }
  

  
  /*allocate additional vectors for computation*/
  *var = (double*) calloc(sizeof(double), (num_cells));
  *cgup = (double*) calloc(sizeof(double), (num_cells));
  *cnorm = (double*) calloc(sizeof(double), (num_internal_cells));
}//memory allocation

void write_vtk(char *file_in, char *scalars_name, int *local_global_index, int num_internal_cells, double *scalars, char *part_type, int myrank) 
{
  int i;
  char file_vtk_out_name [100], buf[10];
  strcpy(file_vtk_out_name, file_in);
  // Finding ".geo.bin" and extracting just the needed part
  char * pch;
  pch = strstr (file_vtk_out_name,".geo.bin");
  strncpy (pch,".",8);
  strcat(file_vtk_out_name, scalars_name);
  strcat(file_vtk_out_name, ".");
  strcat(file_vtk_out_name, part_type);
  strcat(file_vtk_out_name, ".rank");
  snprintf (buf, 10, "%d.vtk", myrank);
  strcat(file_vtk_out_name, buf);
  printf("\n%s\n", file_vtk_out_name);
  
  test_distribution(file_in, file_vtk_out_name, local_global_index, num_internal_cells, scalars);
}//write_vtk


