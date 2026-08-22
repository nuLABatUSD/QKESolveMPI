#!/usr/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: bash test_efficiencyRK.sh N_cores"
    echo "example:"
    echo "    bash test_efficiencyRK.sh 12"
    exit 1
fi

. ./script/script_vars.sh

bash solvingProblems/DistributeWorkLoad.sh 201 20. $1

rm -f coll

mpic++ \
    ${run_code_folder}/run_QKEMPI.cc ${QKE_code} ${MPI_code} \
    -std=c++17 -o coll

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

execute_file="execute_testRK_${1}.sh"

echo "#!/bin/bash" > $execute_file
echo "# FILENAME: 1MeV" >> $execute_file
echo "#SBATCH -A phy240216" >> $execute_file
echo "#SBATCH --nodes=1" >> $execute_file
echo "#SBATCH --ntasks=$1" >> $execute_file
echo "#SBATCH -J RK$1" >> $execute_file
echo "#SBATCH -p shared" >> $execute_file
echo "#SBATCH --time=24:00:00" >> $execute_file
echo "#SBATCH --mail-user=ckishimoto@sandiego.edu" >> $execute_file
echo "#SBATCH --mail-type=all" >> $execute_file

echo "module purge" >> $execute_file
echo "module load modtree/cpu" >> $execute_file
echo "mpiexec -n $1 coll test0 0" >> $execute_file
echo "mpiexec -n $1 coll test1 1" >> $execute_file


echo "Ready to run: sbatch ${execute_file}"

