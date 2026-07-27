#!/usr/bin/bash
if [ $# -ne 3 ]; then
    if [ $# -ne 4 ]; then
        echo "usage: bash compile_QKEMPI-contined.sh continue_file N_steps dN <opt: time, default 12>"
        exit 1
    fi
fi

input_file=$1

if [ ! -f $input_file ]; then
    if [ -f "results/${input_file}" ]; then
        input_file="results/${input_file}"
    fi
    if [ -f "${input_file}.csv" ]; then
        input_file="${input_file}.csv"
    fi
    if [ -f "results/${input_file}.csv" ]; then
        input_file="results/${input_file}.csv"
    fi
fi

if [ ! -f $input_file ]; then
    echo "Error: File ${input_file} not found."
    exit 1
fi

param_file="${input_file%.*}_params.hh"
if [ ! -f $param_file ]; then
    echo "Error: Parameter file ${param_file} does not exist."
    exit 1
fi

if [ $# -eq 4 ]; then
    hrs=$4
else
    hrs="12"
fi

output_base=${input_file%.*}
output_base=${output_base#*/}

. ./script/script_vars.sh 

output_file=$(python3 ${script_folder}/QKEMPI_compile.py $output_base )

numprocs="128"
execute_file="execute_${output_file}.sh"

program_name="${output_file}_run"

if [ -f $program_name ]; then
    rm $program_name
fi

if [ -f $execute_file ]; then
    rm $execute_file
fi

bash solvingProblems/lastlineScript.sh $input_file

if [ ! -f "solvingProblems/restart.hh" ]; then
    echo "lastlineScript.sh didn't work"
    exit 1
fi

mv "solvingProblems/restart.hh" "restart.hh"

cp $param_file run_params.hh

mpic++ ${run_code_folder}/run_QKEMPI_cont.cc ${MPI_code} ${QKE_code} -std=c++11 -o $program_name

if [ -f $program_name ]; then
    echo "#!/bin/bash" > $execute_file
    echo "# FILENAME: 1MeV" >> $execute_file
    echo "#SBATCH -A phy240216" >> $execute_file
    echo "#SBATCH --nodes=1" >> $execute_file
    echo "#SBATCH --ntasks=128" >> $execute_file
    echo "#SBATCH -J ${output_base}" >> $execute_file
    echo "#SBATCH -p wholenode" >> $execute_file
    echo "#SBATCH --time=${hrs}:00:00" >> $execute_file
    echo "#SBATCH --mail-user=ckishimoto@sandiego.edu" >> $execute_file
    echo "#SBATCH --mail-type=all" >> $execute_file
    
    echo "module purge" >> $execute_file
    echo "module load modtree/cpu" >> $execute_file
    echo "mpiexec -n $numprocs $program_name results/${output_file} $2 $3" >> $execute_file

    param_file="results/${output_file}_params.hh"    
    echo $param_file
    git_version=$(git rev-parse --short HEAD )
    cp run_params.hh $param_file
    echo "" >> $param_file
    echo "// Git version ${git_version}" >> $param_file

    echo "Restart of $input_file. Output $2 steps with dN = $3"
    echo "Ready to run for ${hrs} hours: sbatch ${execute_file}"
fi

exit 0



