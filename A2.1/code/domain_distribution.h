/**
 * Domain distribution
 *
 * @date 10-Nov-2014
 * @author M. Bujny, Y. Dong
 */

#ifndef DOMAIN_DISTRIBUTION_H_
#define DOMAIN_DISTRIBUTION_H_

#define MIN(a,b) ((a) < (b) ? (a) : (b))

void allread_calc_global_idx(int** local_global_index, int *nintci_loc, int *nintcf_loc, int *nextci_loc,
                   int *nextcf_loc, int type,  int dual, int nprocs, int myrank,
                   int nintci, int nintcf, int nextci,
                   int nextcf, int** lcc, int* elems, int points_count);
/*
void calc_array_size(int type, int nprocs, int myrank,
                   int nintci, int nintcf, int nextci,
                   int nextcf, int** lcc);
*/
#endif /* DOMAIN_DISTRIBUTION_H_ */