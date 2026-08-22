#!/usr/bin/bash

if [ $# -ne 3 ]; then
    echo "usage: bash run_QKEMPI.sh output_filename"
    exit 1
fi

. ./script/script_vars.sh 

rm -f coll

mpic++ ${run_code_folder}/run_QKEMPI.cc ${MPI_code} ${QKE_code} -std=c++17 -o coll

#mpiexec -n $(num_cores) coll $(outputName) $(newOrOld_collisions)

mpiexec -n $1 coll $2 $3