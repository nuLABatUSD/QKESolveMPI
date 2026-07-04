#!/usr/bin/bash

. ./script/script_vars.sh 

if [ $# -gt 1 ]; then
    echo "usage: bash test_RK_interpolation.sh <optional: test input filename>"
    exit 1
fi

if [ $# -eq 1 ]; then
    if [ -f "test_extrapolation.hh" ]; then
        mv test_extrapolation.hh test_extrapolation-old.hh
    fi
    cp $1 test_extrapolation.hh
fi

output_file="RKextrap"

program_name="${output_file}_run"

if [ -f $program_name ]; then
    rm $program_name
fi

g++ ${run_code_folder}/run_RK_interp.cc ${QKE_code} -std=c++11 -o $program_name

if [ ! -f $program_name ]; then
    echo "g++ error"
    exit 1
fi

if [ -f "${output_file}.csv" ]; then
    echo "mv ${output_file}.csv ${output_file}-old.csv"
    mv ${output_file}.csv ${output_file}-old.csv
fi

./$program_name ${output_file}.csv ${output_file}-eps.csv ${output_file}_original.csv ${output_file}_original-eps.csv
echo "Output printed to ${output_file}.csv and ${output_file}_original.csv"

rm $program_name