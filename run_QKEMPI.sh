#!/usr/bin/bash

rm coll

mpic++ run_QKEMPI.cc QKEMPI.cc collisionsQKE_MPI.cc collisionsQKE.cc QKESolve.cc thermodynamics.cc arrays.cc base_arrays.cc density.cc matrices.cc -std=c++11 -o coll

mpiexec -n 10 coll "results/06-QKEtest"