/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include "util_read_files.h"
#include <string.h>
#include "initialization.h"
int memoryallocation(int ***LCC_local, double **bs_local, double **be_local, double **bn_local, double **bw_local, double **bh_local, double **bl_local, double **bp_local, double **su_local,int num_internal_cells, int num_cells, int *nintcf, int points_count, double **var, double **cgup, double **oc, double **cnorm);


int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index) {
    /********** START INITIALIZATION **********/
    int i = 0;
    int j = 0;


    /***************read-in the input file*********************/
    /*for oneread only processor 0 reads the data*/
    if (strcmp(read_type, "oneread") == 0 ) {
        if(0 == myrank){
        int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                       &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                                       &*points, &*elems);
	if ( f_status != 0 ) return f_status;
        }                               
                                       }
    /*for allread all processors read the data*/
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
   int /*sendcount, recvcount,*/ source;
   int Nintci_loc, Nintcf_loc, Nextci_loc, Nextcf_loc;/*local index for the beginning and ending of internal/external cells*/

/*initialize the array which will only be contained locally*/ 
      int **LCC_local; 
      double *bs_local, *be_local, *bn_local, *bw_local, *bh_local, *bl_local;
      double *bp_local; 
      double *su_local; 


   /*int length_loc_index;*/
  /*int local_global_index_temp*/
   
        /*sendcount = 1;
        recvcount = 1;
        source = 0;
       MPI_Scatter(&nintci_loc[0],sendcount,MPI_INT,&Nintci_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);
       MPI_Scatter(&nintcf_loc[0],sendcount,MPI_INT,&Nintcf_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);
       MPI_Scatter(&nextci_loc[0],sendcount,MPI_INT,&Nextci_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);
       MPI_Scatter(&nextcf_loc[0],sendcount,MPI_INT,&Nextcf_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);*/
       
   
   /*if(strcmp(read_type, "oneread") == 0){
      
      MPI_Status Status[6];
      
      if (0 == rank){
      allread_calc_global_idx( &*local_global_index, &Nintci_loc, &Nintcf_loc, &Nextci_loc,
			     int *nextcf_loc, int type, int dual, int nprocs, int myrank,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count);
          for (int dest = 1; dest<nprocs; dest++){
            
            int XXXXXXXXXXXXXXXXXXX A BUNCH OF TEMP something For passing the data
            allread_calc_global_idx(int** local_global_index, int *nintci_loc, int *nintcf_loc, int *nextci_loc,
			     int *nextcf_loc, int type, int dual, int nprocs, int myrank,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count);
            MPI_Send(&length_loc_index_temp,1,MPI_INT,dest,dest,MPI_COMM_WORLD);
            MPI_Send(&Nintci_loc_temp,1,MPI_INT,dest,dest,MPI_COMM_WORLD);
            MPI_Send(&Nintcf_loc_temp,1,MPI_INT,dest,dest,MPI_COMM_WORLD);
            MPI_Send(&Nextci_loc_temp,1,MPI_INT,dest,dest,MPI_COMM_WORLD);
            MPI_Send(&Nextcf_loc_temp,1,MPI_INT,dest,dest,MPI_COMM_WORLD);
            MPI_Send(&local_global_index_temp[0],length_loc_index_temp,MPI_INT,dest,dest,MPI_COMM_WORLD);
            
          }
      }
      
      if (rank>0){
      MPI_Recv(&length_loc_index,1,MPI_INT,0,rank,MPI_COMM_WORLD,Status[4]);
      if ((local_global_index = (int*)malloc(length_loc_index*sizeof(int)))== NULL){
      fprintf(stderr, "malloc(local_global_index) failed\n");
        return -1;
      }
      MPI_Recv(&Nintci_loc,1,MPI_INT,0,rank,MPI_COMM_WORLD,Status[0]);
      MPI_Recv(&Nintcf_loc,1,MPI_INT,0,rank,MPI_COMM_WORLD,Status[1]);
      MPI_Recv(&Nextci_loc,1,MPI_INT,0,rank,MPI_COMM_WORLD,Status[2]);
      MPI_Recv(&Nextcf_loc,1,MPI_INT,0,rank,MPI_COMM_WORLD,Status[3]);
      MPI_Recv(&local_global_index[0],length_loc_index,MPI_INT,0,rank,MPI_COMM_WORLD,Status[5]);
      
       
      
      
      }
   }// strcmp(read_type, "oneread") == 0*/
   

   if(strcmp(read_type, "allread") == 0){/**** 1****/
     
allread_calc_global_idx( &(*local_global_index), &Nintci_loc, &Nintcf_loc, &Nextci_loc, &Nextcf_loc, part_type, read_type, nprocs, myrank,/*2,0,*/ *nintci, *nintcf, *nextci, *nextcf, *lcc, *elems, *points_count);     

      num_cells = Nextcf_loc - Nintci_loc +1;
      num_internal_cells = Nintcf_loc - Nintci_loc +1;

    /************************ Array memory allocation *******************************/
        
    memoryallocation(&LCC_local, &bs_local, &be_local, &bn_local, &bw_local, &bh_local, &bl_local, &bp_local, &su_local, /**points_count, &points_local, &elems_local,*/ num_internal_cells, num_cells, &*nintcf, *points_count, &*var, &*cgup, &*oc, &*cnorm);


printf("the inside value for nintcf from processor %d is: %d!" ,myrank, num_internal_cells);  

    /*****************   Read Data   ***************/
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
   
}/******  1  ********/

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

/*if(0 == myrank){
printf("LCC_local[17][2] = %d",(LCC_local)[17][2] );

}*/
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

//printf("OK for free no longer used global data lcc from processor %d !\n", myrank ); 
   /* free(su_local_temp);
    free(bp_local_temp);
    free(bh_local_temp);
    free(bl_local_temp);
    free(bw_local_temp);
    free(bn_local_temp);
    free(be_local_temp);
    free(bs_local_temp);*/
   


/*return back element information also to gloabl data */

*nintci = Nintci_loc;
*nintcf = Nintcf_loc;
*nextci = Nextci_loc;
*nextcf = Nextcf_loc;

 


    return 0;
}

int memoryallocation(int ***LCC_local, double **bs_local, double **be_local, double **bn_local, double **bw_local, double **bh_local, double **bl_local, double **bp_local, double **su_local, /*int points_count_local,int ***points_local, int **elems_local,*/ int num_internal_cells, int num_cells, int *nintcf, int points_count, double **var, double **cgup, double **oc, double **cnorm){
    // allocating LCC_local
  int i =0;  
  if ( (*LCC_local = (int**) malloc((num_internal_cells) * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }

    for ( i = 0; i < num_internal_cells + 1; i++ ) {
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
    
 /*
    
    if ( (*elems_local = (int*) malloc((*nintcf + 1) * 8 * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate elems");
        return -1;
    }
    
   
    if ( (*points_local = (int **) calloc(points_count, sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }

    for ( i = 0; i < points_count_local; i++ ) {
        if ( (*points_local[i] = (int *) calloc(3, sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    } */
    
	/*allocate additional vectors for computation*/
    *var = (double*) calloc(sizeof(double), (num_cells));
    *cgup = (double*) calloc(sizeof(double), (num_cells));
    *cnorm = (double*) calloc(sizeof(double), (num_internal_cells));

    
    }//memory allocation

