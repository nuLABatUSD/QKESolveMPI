#!/bin/sh

#. ./script/test_RK_interpolation.sh

#set -e

# =========================
# Change these paths
# =========================

#RK_OUTPUT="RKoutput (1).csv"
read -p "enter RKoutput file name: " RK_OUTPUT
#this file must be inside the folder solvingProblems


#EPSILON_WEIGHTS="06-nu_e_coll-eps.csv"
read -p "enter EPSILON_WEIGHTS file name: " EPSILON_WEIGHTS
#this file must be inside the folder solvingProblems


#EXTRAP_HEADER="test_extrapolation.hh"

#RK_EXTRAP_ORIGINAL="RKextrap_original.csv"
read -p "enter RK_EXTRAP_ORIGINAL file name: " RK_EXTRAP_ORIGINAL
#currently located in the local directory
#aka inside the QKEsolveMPI-Main folder 

#RK_EXTRAP="RKextrap.csv"
read -p "enter RK_EXTRAP file name: " RK_EXTRAP
#currently located in the local directory
#aka inside the QKEsolveMPI-Main folder 

# =========================
# Compile C++ files
#change the name of any files, if the names change
# =========================

echo "Compiling extrapolationCreation..."
g++ ./solvingProblems/extrapolationCreation.cpp -std=c++17 -o extrapolationCreation.cpp


# =========================
# Run pipeline
# =========================

echo "Creating Test_extrapolation.hh..."
./extrapolationCreation.cpp "$RK_OUTPUT"
# will take the rk output csv file to produce test_extrapolation.hh
echo "extrapolationCreation.cpp finished"

echo ""
echo "if you see errors like, error on line 23, 35, and/or 41, do not worry."
echo "the code is still running fine, however it just means it could not find a specific folder"
echo "thus outputting data on the local directory"
echo ""
echo ""

echo "Running RK interpolation script..." 
"./script/test_RK_interpolation.sh" 
# will take in or at least find the test_extrapolation.hh and 
#   run its own code to produce RKextrap_original.csv 
#   and RKextrap.csv
# un-comment the EXTRAP_HEADER incase we need to change file name

echo "RK interpolation finished"

echo "Creating normal RK output graphs..."
python ./solvingProblems/graphingTime.py "$RK_OUTPUT" "$EPSILON_WEIGHTS"
#will produce a bunch of graphs in relating to the RK output
#   and store the files in graphs_output folder

echo "Creating extrapolation comparison graphs..."
python ./solvingProblems/graphingTimeExtrapolation.py "$RK_EXTRAP" "$RK_EXTRAP_ORIGINAL"
#will produce a bunch of graphs relating to the modified version 
#   and the original version of the RKextrap and make a bunch of graphs
#   in the graphs_find_error folder

echo "Done. Graphs created in:"
echo "  graphs_output"
echo "  graphs_find_errors"
