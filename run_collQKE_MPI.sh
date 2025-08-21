#!/usr/bin/bash

rm coll

mpic++ run_collisionsQKE_MPI.cc collisionsQKE_MPI.cc collisionsQKE.cc thermodynamics.cc arrays.cc base_arrays.cc density.cc matrices.cc -std=c++11 -o coll

mpiexec -n 10 coll