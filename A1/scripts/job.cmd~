#!/bin/bash
##
## optional: energy policy tags
##
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = test
#@ node = 1
### schedule the job to 2 to 4 islands 
#@ island_count=1, 1
#@ total_tasks= 1
## other example
##@ tasks_per_node = 16
#@ wall_clock_limit = 0:30:00
##                    0 h 30 min 0 secs
#@ job_name = fire
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/A1/code/
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=mariusz.bujny@tum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
module unload mpi.ibm
module load mpi.intel
module load papi
perf_off

## Running 3 times cojack with text input
mpiexec -n 1 ./gccg text cojack.dat cojack_text_1
mv summary.out summary_cojack_text_1
mpiexec -n 1 ./gccg text cojack.dat cojack_text_2
mv summary.out summary_cojack_text_2
mpiexec -n 1 ./gccg text cojack.dat cojack_text_3
mv summary.out summary_cojack_text_3
## Running 3 times tjunc with text input
mpiexec -n 1 ./gccg text tjunc.dat tjunc_text_1
mv summary.out summary_tjunc_text_1
mpiexec -n 1 ./gccg text tjunc.dat tjunc_text_2
mv summary.out summary_tjunc_text_2
mpiexec -n 1 ./gccg text tjunc.dat tjunc_text_3
mv summary.out summary_tjunc_text_3

## Running 3 times cojack with binary input
mpiexec -n 1 ./gccg bin cojack.bin cojack_bin_1
mv summary.out summary_cojack_bin_1
mpiexec -n 1 ./gccg bin cojack.bin cojack_bin_2
mv summary.out summary_cojack_bin_2
mpiexec -n 1 ./gccg bin cojack.bin cojack_bin_3
mv summary.out summary_cojack_bin_3
## Running 3 times tjunc with binary input
mpiexec -n 1 ./gccg bin tjunc.dat tjunc_bin_1
mv summary.out summary_tjunc_bin_1
mpiexec -n 1 ./gccg bin tjunc.dat tjunc_bin_2
mv summary.out summary_tjunc_bin_2
mpiexec -n 1 ./gccg bin tjunc.dat tjunc_bin_3
mv summary.out summary_tjunc_bin_3