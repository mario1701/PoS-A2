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

export SCOREP_ENABLE_TRACING=false
export SCOREP_ENABLE_PROFILING=true

## Running 3 times cojack with text input
#mpiexec -n 8 ./gccg pent.geo.bin dual allread
#mv scorep* profiling_8_pent_dual_allread

mpiexec -n 8 ./gccg drall.geo.bin dual allread
mv scorep* profiling_8_drall_dual_allread

mpiexec -n 8 ./gccg drall.geo.bin dual oneread
mv scorep* profiling_8_drall_dual_oneread

mpiexec -n 8 ./gccg cojack.geo.bin dual allread
mv scorep* profiling_8_cojack_dual_allread

mpiexec -n 8 ./gccg cojack.geo.bin dual oneread
mv scorep* profiling_8_cojack_dual_oneread

export SCOREP_ENABLE_TRACING=true
export SCOREP_ENABLE_PROFILING=false
export SCOREP_METRIC_PAPI=PAPI_DP_OPS
#module load papi


## Running 3 times cojack with text input
#mpiexec -n 8 ./gccg pent.geo.bin dual allread

mpiexec -n 8 ./gccg drall.geo.bin dual allread
mv scorep* tracing_8_drall_dual_allread

mpiexec -n 8 ./gccg drall.geo.bin dual oneread
mv scorep* tracing_8_drall_dual_oneread

mpiexec -n 8 ./gccg cojack.geo.bin dual allread
mv scorep* tracing_8_cojack_dual_allread

mpiexec -n 8 ./gccg cojack.geo.bin dual oneread
mv scorep* tracing_8_cojack_dual_oneread
