/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdlib.h>

#include "util_read_files.h"
#include "initialization.h"

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
    if (*read_type == "oneread") {
        if(0 == myrank){
        int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                       &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                                       &*points, &*elems);
        }                               
                                       }
    /*for allread all processors read the data*/
    if (*read_type == "allread"){
        int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                       &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                                       &*points, &*elems);
    }

    if ( f_status != 0 ) return f_status;

   /*********************decide how to split the data*************************/
   
   
   
   
   
   
   
   
   
   
   /************************data distribution*************************/
   
   
   int num_cells;
   int num_internal_cells, num_internal_points;
   int sendcount, recvcount, source;
   int Nintci_loc, Nintcf_loc, Nextci_loc, Nextcf_loc;
        /*sendcount = 1;
        recvcount = 1;
        source = 0;
       MPI_Scatter(&nintci_loc[0],sendcount,MPI_INT,&Nintci_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);
       MPI_Scatter(&nintcf_loc[0],sendcount,MPI_INT,&Nintcf_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);
       MPI_Scatter(&nextci_loc[0],sendcount,MPI_INT,&Nextci_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);
       MPI_Scatter(&nextcf_loc[0],sendcount,MPI_INT,&Nextcf_loc,recvcount,MPI_INT,source,MPI_COMM_WORLD);*/
       
   
   /*for allread*/
   
      allread_calc_global_idx(................);
      
      /*suppose I get an index array called local_global_index*/ 
      num_cells = Nextcf_loc - Nintci_loc +1;
      num_internal_cells = Nintcf_loc - Nintci_loc +1;
      num_points = 8 * num_internal_cells;
 
 
      int **LCC_local; 
      double *bs_local, *be_local, *bn_local, *bw_local, *bh_local, *bl_local;
      double *bp_local; 
      double *su_local;   
      
      int points_count_local;    
      int** points_local;    /// coordinates of the points that define the cells - size [points_cnt][3]
      int* elems_local;    /// definition of the cells using their nodes (points) - each cell has 8 points

    /************************ Array memory allocation *******************************/
        
    // allocating LCC_local
    if ( (LCC_local = (int**) malloc((num_internal_cells) * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }

    for ( i = 0; i < num_internal_cells + 1; i++ ) {
        if ( ((LCC_local)[i] = (int*) malloc(6 * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
            return -1;
        }
    }
    
    /*read LCC for LCC_local*/
    for (i =Nintci_loc; i <=Nintcf_loc; i++){
        for (j = 0; j<6; j++){
        LCC_local[i][j] = *LCC[local_global_index[i]][j];
        }
    }
    
    // allocate other arrays
    if ( (bs_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BS) failed\n");
        return -1;
    }

    if ( (be_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BE) failed\n");
        return -1;
    }

    if ( (bn_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BN) failed\n");
        return -1;
    }

    if ( (bw_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BW) failed\n");
        return -1;
    }

    if ( (bh_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BL) failed\n");
        return -1;
    }

    if ( (bl_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BH) failed\n");
        return -1;
    }

    if ( (bp_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BP) failed\n");
        return -1;
    }

    if ( (su_local = (double *) malloc((num_cells) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(SU) failed\n");
        return -1;
    }

    /*read arrays*/
    for (i =Nintci_loc; i <=Nextcf_loc; i++){
        bs_local[i] = *bs[local_global_index[i]];
        be_local[i] = *be[local_global_index[i]];
        bn_local[i] = *bn[local_global_index[i]];
        bw_local[i] = *bw[local_global_index[i]];
        bh_local[i] = *bh[local_global_index[i]];
        bl_local[i] = *bl[local_global_index[i]];
        bp_local[i] = *bp[local_global_index[i]];
        su_local[i] = *su[local_global_index[i]];
    }
   
   /************************ read geometry   ***********************************/
    // allocate elems
    
    if ( (elems_local = (int*) malloc(num_internal_cells * 8 * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate elems");
        return -1;
    }

    // read elems
    for ( i = Nintci_loc; i < (num_internal_cells * 8); i++ ) {
        elems_local[];
    }

    fread(points_count, sizeof(int), 1, fp);

   
   
   
   
   
   /*************need to consider about how many memory do we need!!!**************/
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
