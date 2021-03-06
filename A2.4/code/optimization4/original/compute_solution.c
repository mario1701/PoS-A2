/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */

//#define OVERLAPPING

//#define SCOREP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#ifdef SCOREP
#include <scorep/SCOREP_User.h>
#endif

int compute_solution(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int nghb_cnt, 
                     int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst, int int_access_cnt, int ext_access_cnt, int* int_access, int* ext_access){
  
  // Add SCOREP manual instrumentation
  #ifdef SCOREP
  SCOREP_USER_REGION_DEFINE(handle1);
  SCOREP_USER_REGION_DEFINE(handle2);
  SCOREP_USER_REGION_DEFINE(handle3);
  SCOREP_USER_REGION_DEFINE(handle4);
  SCOREP_USER_REGION_DEFINE(handle5);
  SCOREP_USER_REGION_DEFINE(handle6);
  SCOREP_USER_REGION_DEFINE(handle7);
  SCOREP_USER_REGION_DEFINE(handle8);
  SCOREP_USER_REGION_DEFINE(handle9);
  SCOREP_USER_REGION_DEFINE(handle10);
  SCOREP_USER_REGION_DEFINE(handle_break);
  #endif
  
  #ifdef SCOREP
  SCOREP_USER_REGION_BEGIN( handle1, "handle1 - Initialization of variables and reference residuals.",SCOREP_USER_REGION_TYPE_COMMON );
  #endif
  
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;
    
    /** the reference residual */
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }
    
    #ifdef SCOREP
    SCOREP_USER_REGION_END( handle1 );
    SCOREP_USER_REGION_BEGIN( handle2, "handle2 - 1st Allreduce.",SCOREP_USER_REGION_TYPE_COMMON );
    #endif
    
    // A2.3
    double global_resref = 0;
    MPI_Allreduce(&resref, &global_resref, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    resref = global_resref;
    
    #ifdef SCOREP
    SCOREP_USER_REGION_END( handle2 );
    SCOREP_USER_REGION_BEGIN( handle3, "handle3 - Calculation of the residue sum.",SCOREP_USER_REGION_TYPE_COMMON );
    #endif
    

    resref = sqrt(resref);
    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }
    
    #ifdef SCOREP
    SCOREP_USER_REGION_END( handle3 );
    SCOREP_USER_REGION_BEGIN( handle4, "handle4 - Memory allocation.",SCOREP_USER_REGION_TYPE_COMMON );
    #endif

    
    // Counting the number of ghost cells to extend the direc1
    int ghost_cells_recv = 0, ghost_cells_send = 0;
    int proc, i, j;
    
    for (proc = 0; proc < nghb_cnt; proc++) {
      ghost_cells_recv += recv_cnt[proc];
      ghost_cells_send += send_cnt[proc];
    }
    
    
    /** the computation vectors */
    // TODO:
    double *direc1 = (double *) calloc(sizeof(double), ((nextcf + 1) + ghost_cells_recv));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    
    // Determine displacements for sending
    int **displacements = (int **) malloc(sizeof(double)*nghb_cnt);
    int **blocklenghts = (int **) malloc(sizeof(double)*nghb_cnt);
    
    for(proc = 0; proc < nghb_cnt; proc++) {
      displacements[proc] = (int*)calloc(send_cnt[proc],sizeof(int));
      blocklenghts[proc] = (int*)calloc(send_cnt[proc],sizeof(int));
    }
    
    j = 0;
    for (proc = 0; proc < nghb_cnt; proc++) {
      for (i = 0; i < send_cnt[proc]; i++) {
	displacements[proc][i] = global_local_index[send_lst[proc][i]];
	blocklenghts[proc][i] = 1;
      }
    }
    
    MPI_Request *request1, *request2;
    request1 = (MPI_Request *) malloc(sizeof(*request1)*nghb_cnt);
    request2 = (MPI_Request *) malloc(sizeof(*request2)*nghb_cnt);
    
    MPI_Status status;
    
    MPI_Datatype *indextype;
    indextype = (MPI_Datatype *) malloc(sizeof(*indextype)*nghb_cnt);
    
    for (proc = 0; proc < nghb_cnt; proc++) {
      MPI_Type_indexed(send_cnt[proc], blocklenghts[proc], displacements[proc], MPI_DOUBLE, &(indextype[proc]));
      MPI_Type_commit(&(indextype[proc]));
    }
    
    #ifdef SCOREP
    SCOREP_USER_REGION_END( handle4 );
    SCOREP_USER_REGION_BEGIN( handle5, "handle5 - Computation phase1. direc1 update.",SCOREP_USER_REGION_TYPE_COMMON );
    #endif


    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

	  #ifdef SCOREP
	  SCOREP_USER_REGION_END( handle5 );
	  SCOREP_USER_REGION_BEGIN( handle6, "handle6 - Computation phase1. direc1 communication",SCOREP_USER_REGION_TYPE_COMMON );
	  #endif
        
	  // Communication of direc1 - start
	  
	  for (proc = 0; proc < nghb_cnt; proc++) {
	    MPI_Isend(direc1, 1, indextype[proc], nghb_to_rank[proc], 0, MPI_COMM_WORLD, &(request1[proc]));
	  }
	  
	  if (iter>1) {
	  // prepare additional arrays for the next iteration step
		  if ( nor == nomax ) {
		      nor = 1;
		  } else {
		      if ( nor == 1 ) {
			  for ( nc = nintci; nc <= nintcf; nc++ ) {
			      adxor1[nc] = direc2[nc];
			  }
		      } else {
			  if ( nor == 2 ) {
			      for ( nc = nintci; nc <= nintcf; nc++ ) {
				  adxor2[nc] = direc2[nc];
			      }
			  }
		      }

		      nor++;
		  }
		  nor1 = nor - 1;
	  }
	  
	  // Reference position in the direc1
	  int ref_pos = nextcf + 1;
	  
	  for (proc = 0; proc < nghb_cnt; proc++) {
 
	    #ifdef OVERLAPPING
 	    MPI_Irecv(&(direc1[ref_pos]), recv_cnt[proc], MPI_DOUBLE, nghb_to_rank[proc], 0, MPI_COMM_WORLD, &(request2[proc]));
	    #else
	    MPI_Recv(&(direc1[ref_pos]), recv_cnt[proc], MPI_DOUBLE, nghb_to_rank[proc], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    #endif
	    
	    ref_pos += recv_cnt[proc];

	  }
	 
	  // Communication of direc1 - stop
	  
	  #ifdef SCOREP
	  SCOREP_USER_REGION_END( handle6 );
	  SCOREP_USER_REGION_BEGIN( handle7, "handle7 - Computation phase1. direc2 computation",SCOREP_USER_REGION_TYPE_COMMON );
	  #endif
        
	#ifdef OVERLAPPING
	for ( i = 0; i < int_access_cnt; i++ ) {
	nc = int_access[i];
	#else
	for ( nc = nintci; nc <= nintcf; nc++ ) {
	#endif

		
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[global_local_index[lcc[nc][0]]]
                         - be[nc] * direc1[global_local_index[lcc[nc][1]]] - bn[nc] * direc1[global_local_index[lcc[nc][2]]]
                         - bw[nc] * direc1[global_local_index[lcc[nc][3]]] - bl[nc] * direc1[global_local_index[lcc[nc][4]]]
                         - bh[nc] * direc1[global_local_index[lcc[nc][5]]];
			
			 
        }

        #ifdef OVERLAPPING
 	  for (proc = 0; proc < nghb_cnt; proc++) {
 		MPI_Wait(&(request2[proc]), &status);
 	  }
	
	for ( i = 0; i < ext_access_cnt; i++ ) {
	nc = ext_access[i];

		
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[global_local_index[lcc[nc][0]]]
                         - be[nc] * direc1[global_local_index[lcc[nc][1]]] - bn[nc] * direc1[global_local_index[lcc[nc][2]]]
                         - bw[nc] * direc1[global_local_index[lcc[nc][3]]] - bl[nc] * direc1[global_local_index[lcc[nc][4]]]
                         - bh[nc] * direc1[global_local_index[lcc[nc][5]]];
			
			 
        }
        #endif
      
        
        /********** END COMP PHASE 1 **********/
	
	  #ifdef SCOREP
	  SCOREP_USER_REGION_END( handle7 );
	  SCOREP_USER_REGION_BEGIN( handle8, "handle8 - Computation phase2. occ computation",SCOREP_USER_REGION_TYPE_COMMON );
	  #endif

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for ( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }
            
            // A2.3
            double global_occ = 0.0;
	    MPI_Allreduce(&occ, &global_occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    occ = global_occ;

            oc1 = occ / cnorm[1];
            for ( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor1[nc];
                }
                
                // A2.3
		double global_occ = 0.0;
		MPI_Allreduce(&occ, &global_occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		occ = global_occ;
		

                oc1 = occ / cnorm[1];
                oc2 = 0;
                occ = 0;
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor2[nc];
                }
                
                // A2.3
	
		MPI_Allreduce(&occ, &global_occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		occ = global_occ;

                oc2 = occ / cnorm[2];
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                }

                if2++;
            }
        }
        
       #ifdef SCOREP
      SCOREP_USER_REGION_END( handle8 );
      SCOREP_USER_REGION_BEGIN( handle9, "handle9 - Computation phase2. residual_ratio computation - before break",SCOREP_USER_REGION_TYPE_COMMON );
      #endif

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }
        
	// A2.3
	double global_cnorm_nor = 0.0, global_omega = 0.0;
	MPI_Allreduce(&(cnorm[nor]), &global_cnorm_nor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&omega, &global_omega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	cnorm[nor] = global_cnorm_nor;
	omega = global_omega;

        omega = omega / cnorm[nor];
        double res_updated = 0.0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }
        
	// A2.3
	double global_res_updated = 0.0;
	MPI_Allreduce(&res_updated, &global_res_updated, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	res_updated = global_res_updated;

        res_updated = sqrt(res_updated);
        *residual_ratio = res_updated / resref;

       #ifdef SCOREP
      SCOREP_USER_REGION_END( handle9 );
      #endif
        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

      #ifdef SCOREP
      SCOREP_USER_REGION_BEGIN( handle_break, "handle9 - Computation phase2. residual_ratio computation - after break",SCOREP_USER_REGION_TYPE_COMMON );
      #endif

      iter++;

      int nor_temp;
      nor_temp = nor;
      
      // prepare additional arrays for the next iteration step
      if ( nor == nomax ) {
	  nor = 1;
      } else {
	  if ( nor == 1 ) {
	      for ( nc = nintci; nc <= nintcf; nc++ ) {
		  dxor1[nc] = direc1[nc];
	      }
	  } else {
	      if ( nor == 2 ) {
		  for ( nc = nintci; nc <= nintcf; nc++ ) {
		      dxor2[nc] = direc1[nc];
		  }
	      }
	  }

	  nor++;
      }

      nor = nor_temp;

        /********** END COMP PHASE 2 **********/
    }
    
      #ifdef SCOREP
      SCOREP_USER_REGION_END( handle_break );
      SCOREP_USER_REGION_BEGIN( handle10, "handle10 - Memory freeing",SCOREP_USER_REGION_TYPE_COMMON );
      #endif

    for (i = 0; i < nghb_cnt; i++){
      free(displacements[i]);
    }

    free(displacements);
    
    free(indextype);
    
    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);
    free(resvec);
    
    free(request1);
    free(request2);

    return iter;
    
      #ifdef SCOREP
      SCOREP_USER_REGION_END( handle10 );
      #endif

    

}


