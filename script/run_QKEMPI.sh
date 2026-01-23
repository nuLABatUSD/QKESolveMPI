#!/usr/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: bash run_QKEMPI.sh output_filename"
    exit 1
fi

. ./script/script_vars.sh 

rm coll

mpic++ ${run_code_folder}/run_QKEMPI.cc ${MPI_code} ${QKE_code} -std=c++11 -o coll

mpiexec -n 10 coll $1