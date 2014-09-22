#!/bin/bash

#@ class            = clallmds+
#@ job_name         = job
#@ total_tasks      = 1
#@ node             = 1
#@ node_usage       = not_shared
#@ wall_clock_limit = 20:00:00
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).log
#@ job_type         = mpich
#@ environment      = COPY_ALL 
#@ queue

module load gnu-env/4.9.0
module load fftw/3.3.3_gnu49

JOB=$(head -1 input/dft.in | cut -c11-13)
./../mdft
cp output/rdf.out rdf.$JOB.out
rm -rf output
rm input/ck*
rm x.sh
