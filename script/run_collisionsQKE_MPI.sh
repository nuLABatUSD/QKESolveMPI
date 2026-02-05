#!/usr/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: bash run_collisionsQKE_MPI.sh N_cores"
    exit 1
fi

. ./script/script_vars.sh 

rm coll

mpic++ ${run_code_folder}/run_collisionsQKE_MPI.cc ${QKE_code} ${MPI_code} -std=c++11 -o coll

mpiexec -n $1 coll