#!/bin/bash

# Configure Extrae
export EXTRAE_HOME=/usr/local/software/jureca/Stages/2016a/software/Extrae/3.3.0-iimpi-8.2.5-GCC-4.9.3-2.25
export EXTRAE_CONFIG_FILE=./extrae.xml

# Load the tracing library (choose C/Fortran)
export LD_PRELOAD=$EXTRAE_HOME/lib/libseqtrace.so
#export LD_PRELOAD=$EXTRAE_HOME/lib/libmpitracef.so

# Run the program
$*
