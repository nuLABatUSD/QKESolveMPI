#!/usr/bin/bash

. ./script/script_vars.sh 

rm coll

mpic++ ${run_code_folder}/test_QKEMPI_Coll.cc ${MPI_code} ${QKE_code} -std=c++11 -o coll

mpiexec -n 10 coll