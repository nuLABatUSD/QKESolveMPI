#!/usr/bin/bash

code_folder="code"
run_code_folder="run"
script_folder="script"

MPI_code="${code_folder}/QKEMPI.cc ${code_folder}/collisionsQKE_MPI.cc"
coherent_code="${code_folder}/QKESolve.cc ${code_folder}/thermodynamics.cc ${code_folder}/arrays.cc ${code_folder}/density.cc ${code_folder}/matrices.cc ../base_code/base_arrays.cc"
QKE_code="${code_folder}/collisionsQKE.cc ${coherent_code}"
