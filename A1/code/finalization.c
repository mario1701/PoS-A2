/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <string.h>
#include "util_write_files.h"
#include "vol2mesh.h"

void finalization(char* file_in, int total_iters, double residual_ratio,
		  int nintci, int nintcf, double* var, double* cgup, double* su, int **lcc, char *prefix, long long *counters, float rtime, float ptime, float mflops, long long flpops, long long *countersm)
{
  
  char file_out[30] = "summary.out";
  char file_su_vtk[30];
  char file_var_vtk[30];
  char file_cgup_vtk[30];
  
  strcpy(file_su_vtk, prefix); /* copy name into the new var */
  strcat(file_su_vtk, ".SU.vtk"); /* add the extension */
  strcpy(file_var_vtk, prefix); /* copy name into the new var */
  strcat(file_var_vtk, ".VAR.vtk"); /* add the extension */
  strcpy(file_cgup_vtk, prefix); /* copy name into the new var */
  strcat(file_cgup_vtk, ".CGUP.vtk"); /* add the extension */
  
  int status = write_result(file_in, file_out, nintci, nintcf, var, total_iters, residual_ratio, counters, rtime, ptime, mflops, flpops, countersm);
  
  if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);
  
  // Writing VTK files
  int nodeCnt;
  int **points, **elems;
  
  vol2mesh(nintci, nintcf, lcc, &nodeCnt, &points, &elems);
  write_result_vtk(file_su_vtk, nintci, nintcf, nodeCnt, points, elems, su);
  write_result_vtk(file_var_vtk, nintci, nintcf, nodeCnt, points, elems, var);
  write_result_vtk(file_cgup_vtk, nintci, nintcf, nodeCnt, points, elems, cgup);

  
}
		  
		  
