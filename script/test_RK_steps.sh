#!/usr/bin/bash

. ./script/script_vars.sh 

if [ $# -gt 1 ]; then
    echo "usage: bash test_RK_steps.sh <optional: test cases filename>"
    exit 1
fi

if [ $# -eq 1 ]; then
    if [ -f "test_cases.hh" ]; then
        mv test_cases.hh test_cases-old.hh
    fi
    cp $1 test_cases.hh
fi

output_file="RKtest"

program_name="${output_file}_run"

if [ -f $program_name ]; then
    rm $program_name
fi

mpic++ ${run_code_folder}/run_RK_test.cc ${MPI_code} ${QKE_code} -std=c++11 -o $program_name

if [ ! -f $program_name ]; then
    echo "mpic++ error"
    exit 1
fi

if [ -f "RKoutput.csv" ]; then
    echo "mv RKoutput.csv RKoutput-old.csv"
    mv RKoutput.csv RKoutput-old.csv
fi

mpiexec -n 10 $program_name RKoutput.csv
echo "Output printed to RKoutput.csv"

rm $program_name