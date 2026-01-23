#!/usr/bin/bash

. ./script/script_vars.sh 

rm coll

mpic++ ${run_code_folder}/run_collisionsQKE_MPI.cc ${QKE_code} ${MPI_code} -std=c++11 -o coll

mpiexec -n 10 coll