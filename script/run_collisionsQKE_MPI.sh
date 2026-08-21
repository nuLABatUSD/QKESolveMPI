#!/usr/bin/bash

if [ $# -ne 2 ]; then
    if [ $# -ne 3]; then
        echo "usage: bash run_collisionsQKE_MPI.sh N_cores <repeat default=1>"
        echo "example:"
        echo "    bash run_collisionsQKE_MPI.sh 12 2"
        exit 1
    fi
fi

if [ $# -eq 2 ]; then
    m="$1 ./coll"
fi
if [ $# -eq 3 ]; then
    m="$1 ./coll $2"
fi

N_CORES="$1"
MODE="$2"

. ./script/script_vars.sh


#echo
#echo "Select collision test:"
#echo "    0 = original collision constructor"
#echo "    1 = optimized collision constructor"
#echo "    2 = run both and compare"
#echo


#read -p "Enter mode [0/1/2]: " MODE

#case "$MODE" in
#    0|1|2)
#        ;;
#    *)
#        echo "Error: mode must be 0, 1, or 2."
#        exit 1
#        ;;
#esac

#this is here incase we need to select 

#MODE="1"

rm -f coll

mpic++ \
    ${run_code_folder}/run_collisionsQKE_MPI.cc ${QKE_code} ${MPI_code} \
    -std=c++17 -o coll

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

mpiexec -n $m
