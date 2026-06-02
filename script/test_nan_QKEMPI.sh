#!/usr/bin/bash
if [ $# -ne 2 ]; then
    echo "usage: bash test_nan_QKEMPI.sh output_filename input_filename"
    exit 1
fi

module purge
module load modtree/cpu
module load python/3.9.5

output_base=$1
input_base=$2

. ./script/script_vars.sh 

if [ -f "density_test.hh" ]; then
    rm density_test.hh
fi

python3 ${script_folder}/QKEMPI_nan_test.py $input_base
output_file=$(python3 ${script_folder}/QKEMPI_compile.py $output_base )

if [ ! -f "density_test.hh" ]; then
    echo "QKEMPI_nan_test.py did not create density_test.hh"
    exit 1
fi

cp ${input_base}_params.hh run_params.hh

numprocs="128"
execute_file="execute_${output_file}.sh"

program_name="${output_file}_run"

if [ -f $program_name ]; then
    rm $program_name
fi

if [ -f $execute_file ]; then
    rm $execute_file
fi

mpic++ ${run_code_folder}/run_QKEMPI_test_nan.cc ${MPI_code} ${QKE_code} -std=c++11 -o $program_name

if [ -f $program_name ]; then
    echo "#!/bin/bash" > $execute_file
    echo "# FILENAME: 1MeV" >> $execute_file
    echo "#SBATCH -A phy240216" >> $execute_file
    echo "#SBATCH --nodes=1" >> $execute_file
    echo "#SBATCH --ntasks=128" >> $execute_file
    echo "#SBATCH -J ${output_base}" >> $execute_file
    echo "#SBATCH -p wholenode" >> $execute_file
    echo "#SBATCH --time=01:00:00" >> $execute_file
    echo "#SBATCH --mail-user=ckishimoto@sandiego.edu" >> $execute_file
    echo "#SBATCH --mail-type=all" >> $execute_file
    
    echo "module purge" >> $execute_file
    echo "module load modtree/cpu" >> $execute_file
    echo "mpiexec -n $numprocs $program_name results/${output_file}" >> $execute_file
    echo "rm density_test.hh" >> $execute_file

    param_file="results/${output_file}_params.hh"    
    echo $param_file
    git_version=$(git rev-parse --short HEAD )
    cp run_params.hh $param_file
    echo "" >> $param_file
    echo "// Git version ${git_version}" >> $param_file

    echo "Ready to run: sbatch ${execute_file}"
fi
