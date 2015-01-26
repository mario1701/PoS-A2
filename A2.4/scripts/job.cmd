#!/bin/bash
##
## optional: energy policy tags
##
# DO NOT USE environment = COPY_ALL
#@ job_type = parallel
#@ class = test
#@ node = 1
### schedule the job to 2 to 4 islands 
#@ island_count=1, 1
#@ total_tasks= 8
## other example
##@ tasks_per_node = 8
#@ wall_clock_limit = 0:30:00
##                    0 h 30 min 0 secs
#@ job_name = fire
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/Supercomputers/opt0/
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=never
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
#module unload mpi.ibm
#module load mpi.intel
#perf_off

## Running 3 times cojack with text input
#mpiexec -n 1 ./gccg pent.geo.bin dual allread
#mv scorep* profiling_1_pent_dual_allread
#mpiexec -n 2 ./gccg pent.geo.bin dual allread
#mv scorep* profiling_2_pent_dual_allread
mpiexec -n 8 ./gccg pent.geo.bin dual oneread
#mv scorep* profiling_8_pent_dual_allread
#mpiexec -n 16 ./gccg pent.geo.bin dual allread
#mv scorep* profiling_16_pent_dual_allread
#mpiexec -n 32 ./gccg pent.geo.bin dual allread
#mv scorep* profiling_32_pent_dual_allread


