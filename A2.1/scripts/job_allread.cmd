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
#@ total_tasks= 16
## other example
##@ tasks_per_node = 16
#@ wall_clock_limit = 0:30:00
##                    0 h 30 min 0 secs
#@ job_name = fire
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/Supercomputers/A2/
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

mpiexec -n 4 ./gccg drall.geo.bin classic allread
mv pstats.dat pstats.drall.4.classic.allread.dat

mpiexec -n 9 ./gccg drall.geo.bin classic allread
mv pstats.dat pstats.drall.9.classic.allread.dat

mpiexec -n 12 ./gccg drall.geo.bin classic allread
mv pstats.dat pstats.drall.12.classic.allread.dat

mpiexec -n 4 ./gccg cojack.geo.bin classic allread
mv pstats.dat pstats.cojack.4.classic.allread.dat

mpiexec -n 9 ./gccg cojack.geo.bin classic allread
mv pstats.dat pstats.cojack.9.classic.allread.dat

mpiexec -n 12 ./gccg cojack.geo.bin classic allread
mv pstats.dat pstats.cojack.12.classic.allread.dat

mpiexec -n 4 ./gccg drall.geo.bin dual allread
mv pstats.dat pstats.drall.4.dual.allread.dat

mpiexec -n 9 ./gccg drall.geo.bin dual allread
mv pstats.dat pstats.drall.9.dual.allread.dat

mpiexec -n 12 ./gccg drall.geo.bin dual allread
mv pstats.dat pstats.drall.12.dual.allread.dat

mpiexec -n 4 ./gccg cojack.geo.bin dual allread
mv pstats.dat pstats.cojack.4.dual.allread.dat

mpiexec -n 9 ./gccg cojack.geo.bin dual allread
mv pstats.dat pstats.cojack.9.dual.allread.dat

mpiexec -n 12 ./gccg cojack.geo.bin dual allread
mv pstats.dat pstats.cojack.12.dual.allread.dat

mpiexec -n 4 ./gccg drall.geo.bin nodal allread
mv pstats.dat pstats.drall.4.nodal.allread.dat

mpiexec -n 9 ./gccg drall.geo.bin nodal allread
mv pstats.dat pstats.drall.9.nodal.allread.dat

mpiexec -n 12 ./gccg drall.geo.bin nodal allread
mv pstats.dat pstats.drall.12.nodal.allread.dat

mpiexec -n 4 ./gccg cojack.geo.bin nodal allread
mv pstats.dat pstats.cojack.4.nodal.allread.dat

mpiexec -n 9 ./gccg cojack.geo.bin nodal allread
mv pstats.dat pstats.cojack.9.nodal.allread.dat

mpiexec -n 12 ./gccg cojack.geo.bin nodal allread
mv pstats.dat pstats.cojack.12.nodal.allread.dat
