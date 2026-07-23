#!/bin/sh

csvpath="./solvingProblems/WorkerHeaderCreation"
csvname="$csvpath/load_factors3.csv"

#cal_load_factors
#original file is calculate_diagnostics.cc
diagnosticsCode="./solvingProblems/WorkerHeaderCreation/cal_load_factors.cc"
diagnosticsProgram="./solvingProblems/WorkerHeaderCreation/cal_load_factors.exe"

pathtovars="../base_code-main/base_arrays.cc"
densitypath="./code/density.cc"
arrpath="./code/arrays.cc"
qkepath="./code/collisionsQKE.cc"
matrixpath="./code/matrices.cc"
thermopath="./code/thermodynamics.cc"
#the relative pathways here may need to change 
#   depending on where the base code main files are located 
#read -r -p "Enter collisionType: " collisionType
#read -r -p "Enter N_TRAP name: " N_TRAP
#read -r -p "Enter EPS_MAX file name: " EPS_MAX
#read -r -p "Enter N_CORES file name: " N_CORES

#collisionType="0"
#N_TRAP="201"
#EPS_MAX="20"
#N_CORES="128"

collisionType="$1"
N_TRAP="$2"
EPS_MAX="$3"
N_CORES="$4"

tempCSV="$csvname"
#   this was here, incase i need to see what the printed csv file 
#   looks like and have a reference for later.

if [ ! -d "$csvpath" ]; then
    mkdir -p "$csvpath"
fi
#should create a directory just incase

#==============================================================================
#the code below will compile the files needed to run the diagnostic checks
#which will then create the csv file which we will then for later 

echo
echo "Compiling calculate_diagnostics.cc..."
echo


diagnosticsCode="./solvingProblems/WorkerHeaderCreation/cal_load_factors.cc"
diagnosticsProgram="./solvingProblems/WorkerHeaderCreation/cal_load_factors.exe"
    
g++ \
    "$diagnosticsCode" \
    "$pathtovars" \
    "$densitypath" \
    "$arrpath" \
    "$qkepath" \
    "$matrixpath" \
    "$thermopath" \
    -std=c++17 -O2 -o "$diagnosticsProgram"
#"$compiler" "$diagnosticsSource" -std=c++17 -o "$diagnosticsProgram"


compile_status=$?

echo
echo "Compiler status: $compile_status"

if [ "$compile_status" -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

if [ ! -f "$diagnosticsProgram" ]; then
    echo "Compilation reported success, but this executable was not created:"
    echo "$diagnosticsProgram"
    exit 1
fi

echo
echo "Running diagnostics..."
echo "this may take a few minutes ..." 
echo "Collision type = $collisionType"
echo "N_trap        = $N_TRAP"
echo "eps_max       = $EPS_MAX"
echo "Output CSV    = $tempCSV"
echo


"$diagnosticsProgram" \
    "$collisionType" \
    "$N_TRAP" \
    "$EPS_MAX" \
    "$tempCSV"

run_status=$?

echo
echo "Diagnostics status: $run_status"

if [ "$run_status" -ne 0 ]; then
    echo "calculate_diagnostics.exe failed."
    exit 1
fi

if [ ! -f "$tempCSV" ]; then
    echo "Diagnostics ran, but the CSV was not created:"
    echo "$tempCSV"
    exit 1
fi


echo
echo "Diagnostics completed successfully."
echo "Output saved to:"
echo "$tempCSV"

#=======================================================================

headerSource="./solvingProblems/WorkerHeaderCreation/headerCreation.cc"
headerProgram="./solvingProblems/WorkerHeaderCreation/headerCreation.exe"

echo
echo "Compiling headerCreation.cc..."

g++ \
    "$headerSource" \
    -std=c++17 \
    -O2 \
    -o headerCreation.cc

header_compile_status=$?

echo "Header compiler status: $header_compile_status"

if [ "$header_compile_status" -ne 0 ]; then
    echo "headerCreation.cc failed to compile."
    echo "The CSV will not be deleted."
    #exit 1
fi

if [ ! -f "$headerProgram" ]; then
    echo "The executable was not created:"
    echo "$headerProgram"
    #exit 1
fi


echo
echo "Running header creation..."

./headerCreation.cc \
    "$N_CORES" \
    "$N_TRAP" \
    "$EPS_MAX" \
    "$tempCSV"

header_run_status=$?

echo "Header creation status: $header_run_status"

if [ "$header_run_status" -ne 0 ]; then
    echo "Header creation failed."
    echo "The CSV will not be deleted."
    exit 1
fi

echo
echo "Header file created successfully."

#========================================================================
# 4. Delete the specific file

g++ ./solvingProblems/WorkerHeaderCreation/theRemover.cc -std=c++17 -o theRemover.cc

#this code here deletes the file, this is because the bash command is
#   not working for me
 ./theRemover.cc "$csvname"

echo "CSV file has been successfully deleted."
