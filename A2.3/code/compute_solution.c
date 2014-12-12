/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int compute_solution(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int nghb_cnt, 
                     int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst){
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;

//     printf("Compute solution started\n");
//     int m=0;
//     for (m=0;m<nextcf;m++) {
//      printf("m = %d, nextcf = %d\n", m,  global_local_index[m]);
//     }
    
    /** the reference residual */
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }
    
    // A2.3
    double global_resref = 0;
    MPI_Allreduce(&resref, &global_resref, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    resref = global_resref;
    

    resref = sqrt(resref);
    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }

    
    // Counting the number of ghost cells to extend the direc1
    int ghost_cells_recv = 0, ghost_cells_send = 0;
    int proc, i, j;
    
    for (proc = 0; proc < nghb_cnt; proc++) {
      ghost_cells_recv += recv_cnt[proc];
      ghost_cells_send += send_cnt[proc];
    }
    
    
    /** the computation vectors */
    // TODO: Why calloc the other way around?
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
	//printf("proc = %d, %d\n", nghb_to_rank[proc], displacements[proc][i]);
      }
    }
    
    MPI_Request request;
    MPI_Datatype *indextype;
    indextype = (MPI_Datatype *) malloc(sizeof(*indextype)*nghb_cnt);
    
    for (proc = 0; proc < nghb_cnt; proc++) {
      MPI_Type_indexed(send_cnt[proc], blocklenghts[proc], displacements[proc], MPI_DOUBLE, &(indextype[proc]));
      MPI_Type_commit(&(indextype[proc]));
    }
    
    // Testing
//     for ( nc = nintci; nc <= nintcf; nc++ ) {
// 	printf("nc = %d, %d %d %d %d %d \n", nc, lcc[nc][0], lcc[nc][1], lcc[nc][2], lcc[nc][3], lcc[nc][4], lcc[nc][5], lcc[nc][5]);
//     }


    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

	  // Communication of direc1 - start
	  
	  for (proc = 0; proc < nghb_cnt; proc++) {
// 	    MPI_Type_indexed(send_cnt[proc], blocklenghts[proc], displacements[proc], MPI_DOUBLE, &(indextype[proc]));
// 	    MPI_Type_commit(&(indextype[proc]));
	    //MPI_Send(direc1, 1, indextype, nghb_to_rank[proc], 0, MPI_COMM_WORLD);
	    MPI_Isend(direc1, 1, indextype[proc], nghb_to_rank[proc], 0, MPI_COMM_WORLD, &request);
	  }
	  
	  // Reference position in the direc1
	  int ref_pos = nextcf + 1;
	  
	  for (proc = 0; proc < nghb_cnt; proc++) {
	    MPI_Recv(&(direc1[ref_pos]), recv_cnt[proc], MPI_DOUBLE, nghb_to_rank[proc], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    ref_pos += recv_cnt[proc];
	  }
	  
	  //printf("Comm done\n");
	  
	  // Communication of direc1 - stop
        
        // compute new guess (approximation) for direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
//             direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[nc][0]]
//                          - be[nc] * direc1[lcc[nc][1]] - bn[nc] * direc1[lcc[nc][2]]
//                          - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
//                          - bh[nc] * direc1[lcc[nc][5]];
	  //printf("nextcf = %d\n", (nextcf + 1));
	  //printf("nc = %d, %d %d %d %d %d %d %d \n", nc, global_local_index[lcc[nc][0]], global_local_index[lcc[nc][1]], global_local_index[lcc[nc][2]], global_local_index[lcc[nc][3]], global_local_index[lcc[nc][4]], global_local_index[lcc[nc][5]], lcc[nc][5]);
	  

		
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[global_local_index[lcc[nc][0]]]
                         - be[nc] * direc1[global_local_index[lcc[nc][1]]] - bn[nc] * direc1[global_local_index[lcc[nc][2]]]
                         - bw[nc] * direc1[global_local_index[lcc[nc][3]]] - bl[nc] * direc1[global_local_index[lcc[nc][4]]]
                         - bh[nc] * direc1[global_local_index[lcc[nc][5]]];
			
			 
        }
        
        
        /********** END COMP PHASE 1 **********/

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

        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }

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

    return iter;
}


