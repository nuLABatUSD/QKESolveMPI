#!/usr/bin/bash

# How to use this: input arguments:
# 1: T_cm
# 2-5: k_e, k_mu, k_ebar, k_mubar
# 6: delta m^2
# 7: file stem name: creates files file_run.csv and file_eps.csv
#
# usage: bash run_coherent.sh 32 0.9 1.8 0.9 1.8 2.5e-18 results/coherent_test
# creates files: results/coherent_test_run.csv and results/coherent_test_eps.csv

rm coh

g++ coherentsolve.cc QKESolve.cc thermodynamics.cc arrays.cc base_arrays.cc density.cc -o coh

./coh $1 $2 $3 $4 $5 $6 $7 $8