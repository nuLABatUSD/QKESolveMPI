#!/usr/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: bash test_QKEMPI_Coll.sh N_cores"
    exit 1
fi

. ./script/script_vars.sh 

rm coll

mpic++ ${run_code_folder}/test_QKEMPI_Coll.cc ${MPI_code} ${QKE_code} -std=c++11 -o coll

mpiexec -n $1 coll