#!/usr/bin/bash

if [ $# -ne 4 ]; then
    echo "usage: bash calc_R.sh (0 for nu/nu | 1 for nu_e) N_LS_bins eps_max filename"
    exit 1
fi

. ./script/script_vars.sh 

rm coll

g++ ${run_code_folder}/calculate_diagnostics.cc ${QKE_code} -std=c++11 -o coll

./coll $1 $2 $3 $4