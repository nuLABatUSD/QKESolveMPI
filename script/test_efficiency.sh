#!/usr/bin/bash

if [ $# -ne 1 ]; then
    if [ $# -ne 2 ]; then
        echo "usage: bash run_collisionsQKE_MPI.sh N_cores <repeat default=1>"
        echo "example:"
        echo "    bash run_collisionsQKE_MPI.sh 12 2"
        exit 1
    fi
fi

if [ $# -eq 1 ]; then
    m="$1 ./coll"
fi
if [ $# -eq 2 ]; then
    m="$1 ./coll $2"
fi

. ./script/script_vars.sh

bash solvingProblems/DistributeWorkLoad.sh 201 20. $1

rm -f coll

mpic++ \
    ${run_code_folder}/run_collisionsQKE_MPI.cc ${QKE_code} ${MPI_code} \
    -std=c++17 -o coll

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

execute_file="execute_test_collisions_${1}.sh"

echo "#!/bin/bash" > $execute_file
echo "# FILENAME: 1MeV" >> $execute_file
echo "#SBATCH -A phy240216" >> $execute_file
echo "#SBATCH --nodes=1" >> $execute_file
echo "#SBATCH --ntasks=$1" >> $execute_file
echo "#SBATCH -J ${output_base}" >> $execute_file
echo "#SBATCH -p shared" >> $execute_file
echo "#SBATCH --time=01:00:00" >> $execute_file
echo "#SBATCH --mail-user=ckishimoto@sandiego.edu" >> $execute_file
echo "#SBATCH --mail-type=all" >> $execute_file

echo "module purge" >> $execute_file
echo "module load modtree/cpu" >> $execute_file
echo "mpiexec -n ${m}" >> $execute_file


echo "Ready to run for ${hrs} hours: sbatch ${execute_file}"

