#!/usr/bin/bash

output_base=$1

output_file=$(python3 compile_QKEMPI.py $output_base )

numprocs="128"
execute_file="execute_${output_file}.sh"

program_name="${output_file}_run"

if [ -f $program_name ]; then
    rm $program_name
fi

if [ -f $execute_file ]; then
    rm $execute_file
fi

mpic++ run_QKEMPI.cc QKEMPI.cc collisionsQKE_MPI.cc collisionsQKE.cc QKESolve.cc thermodynamics.cc arrays.cc base_arrays.cc density.cc matrices.cc -std=c++11 -o $program_name

if [ -f $program_name ]; then
    echo "#!/bin/bash" > $execute_file
    echo "# FILENAME: 1MeV" >> $execute_file
    echo "#SBATCH -A phy240216" >> $execute_file
    echo "#SBATCH --nodes=1" >> $execute_file
    echo "#SBATCH --ntasks=128" >> $execute_file
    echo "#SBATCH -J qkesolvef" >> $execute_file
    echo "#SBATCH -p wholenode" >> $execute_file
    echo "#SBATCH --time=12:00:00" >> $execute_file
    echo "#SBATCH --mail-user=ckishimoto@sandiego.edu" >> $execute_file
    echo "#SBATCH --mail-type=all" >> $execute_file
    
    echo "module purge" >> $execute_file
    echo "module load modtree/cpu" >> $execute_file
    echo "mpiexec -n $numprocs $program_name results/${output_file}" >> $execute_file

    cp run_params.hh "results/${output_file}_params.hh"
    
    echo "Ready to run: sbatch ${execute_file}"
fi
