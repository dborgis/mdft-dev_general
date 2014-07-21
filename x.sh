#!/bin/bash

#@ class            = clallmds+
#@ job_name         = dipole
#@ total_tasks      = 1
#@ node             = 1
#@ node_usage       = not_shared
#@ wall_clock_limit = 20:00:00
#@ output           = $(job_name).$(jobid).log
#@ job_type         = mpich
#@ environment      = COPY_ALL 
#@ queue

module load gnu-env/4.7.2
module load openmpi/1.6.3_gnu47
module load fftw/3.3.3_gnu47

JOB=$(head -1 input/dft.in | cut -c11-13)
mpirun ./mdft
cp output/rdf.out rdf.$JOB.out

