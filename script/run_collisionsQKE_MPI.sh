#!/usr/bin/bash
#C:/msys64/usr/bin/bash.exe




#!/usr/bin/bash


if [ $# -ne 2 ]; then
    echo "usage: bash run_collisionsQKE_MPI.sh N_cores distribution_file"
    echo "example:"
    echo "    bash run_collisionsQKE_MPI.sh 12 ../file.hh"
    exit 1
fi


N_CORES="$1"
INPUT_DISTRIBUTION="$2"


DESTINATION="../solvingProblems/generalDistribution.hh"


. ../script/script_vars.sh




# Make sure the input file exists.
if [ ! -f "$INPUT_DISTRIBUTION" ]; then
    echo "Error: distribution file does not exist:"
    echo "    $INPUT_DISTRIBUTION"
    exit 1
fi


# Make sure the destination directory exists.
mkdir -p "../solvingProblems"


# Copy the selected distribution and rename it.
cp "$INPUT_DISTRIBUTION" "$DESTINATION"


if [ $? -ne 0 ]; then
    echo "Error: failed to copy distribution file."
    exit 1
fi


echo "Using distribution:"
echo "    $INPUT_DISTRIBUTION"
echo "Copied to:"
echo "    $DESTINATION"


# Remove the previous executable if it exists.
rm -f coll


# Compile.
mpic++ \
    ${run_code_folder}/run_collisionsQKE_MPI.cc \
    ${QKE_code} \
    ${MPI_code} \
    -std=c++11 \
    -o coll


if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi


# Run.
mpiexec -n "$N_CORES" ./coll




#to run a file on the terminal we can do something like:
#    C:/msys64/usr/bin/bash.exe "c:\Users\ckishimoto\Desktop\QKESolveMPI-extrap\QKESolveMPI-extrap\script\run_collisionsQKE_MPI.sh"
#   this will run the code using bash from mingw and run the directory to the file
#   the PC we are using, or at least this one, has 16 cores
#   simulating the real thing, lets assume 15 cores is what we are dealing with




# C:/msys64/usr/bin/bash.exe "c:\Users\ckishimoto\Desktop\QKESolveMPI-extrap\QKESolveMPI-extrap\script\run_collisionsQKE_MPI.sh"  15
# currently vs code isnt reading the mpi exe and exec file correctly




# this works on command prompt:
#   bash ./run_collisionsQKE_MPI.sh 15


#use for linux:
#   bash run_collisionsQKE_MPI.sh 15






