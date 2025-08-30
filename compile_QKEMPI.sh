#!/usr/bin/bash

set -x
rm coll

numprocs="128"
output_file="QKEthermal.csv"
execute_file="execute_QKEMPI.sh"

sin2theta="0"
deltamsquared="0"

N_step="200"
dN="1"
x_initial="0."
x_final="5.e20"
dx_initial="1.e13"

rm run_params.hh

echo "#define PARAM_DELTA_M_SQUARED $deltamsquared" > run_params.hh
echo "#define PARAM_SIN_2THETA $sin2theta" >> run_params.hh
echo "#define PARAM_N_STEPS $N_step" >> run_params.hh
echo "#define PARAM_DN $dN" >> run_params.hh
echo "#define PARAM_DT_INIT $dx_initial" >> run_params.hh
echo "#define PARAM_T_FINAL $x_final" >> run_params.hh

mpic++ run_QKEMPI.cc QKEMPI.cc collisionsQKE_MPI.cc collisionsQKE.cc QKESolve.cc thermodynamics.cc arrays.cc base_arrays.cc density.cc matrices.cc -std=c++11 -o coll

rm $execute_file

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
echo "mpiexec -n $numprocs coll \"$output_file\" " >> $execute_file
