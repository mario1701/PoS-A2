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
                   int *nextcf_loc, char *part_type, char*read_type, int nprocs, int myrank,
                   int nintci, int nintcf, int nextci,
                   int nextcf, int** lcc, int* elems, int points_count,int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, int **recv_cnt, int*** recv_lst);

void oneread_calc_global_idx(int*** local_global_index, int **nintci_loc, int **nintcf_loc, int **nextci_loc,
			     int **nextcf_loc, char *part_type, char*read_type, int nprocs,
			     int nintci, int nintcf, int nextci,
			     int nextcf, int** lcc, int* elems, int points_count);

#endif /* DOMAIN_DISTRIBUTION_H_ */
