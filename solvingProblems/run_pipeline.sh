#C:/msys64/usr/bin/bash.exe

set -e
pathing="./solvingProblems/"
#change pathing accordingly depending on the locations of where the files are found

#BE WARNED, many things will print onto the terminal when running the script
#   the reason is to see where issues may have occured and to see if the code is
#   comparing things properly

SLURM_OUTPUT="slurm-18475761.txt"
#just the output of all the error cases
ORIGINAL_DATASET="65-fix_10_small.csv"
#this is for the data set with all the density values
MODIFYING_VALUES="06-nu_e_coll-eps.csv"
#this is for the epsilon and weights csv file

# Output files created by programs
ERROR_INFO="error_info.txt"
DATASET_INFO="dataset_info.txt"
EPSILON_WEIGHTS="epsilon_and_weights.txt"
ERROR_FOUND="errorFound.txt"
ULTIMATE_HEADER="UltimateGenericName.hh"
MAXIMUM_HEADER="MaximumGenericName.hh"

#the mode menu:
#cout << "\nWhat do you want to do?\n";
#cout << "1. Find one specific case\n";
#cout << "2. Find multiple specific cases\n";
#cout << "3. Find largest dsm cases\n";
#cout << "Choice: ";
mode="3"


# =========================
# Compile C++ files
# =========================


g++ "${pathing}/Findissues.cpp" -std=c++17 -o Findissues.cpp
g++ "${pathing}/dataCompare.cpp" -std=c++17 -o dataCompare.cpp
g++ "${pathing}/OutputVsData.cpp" -std=c++17 -o OutputVsData.cpp
g++ "${pathing}/testCreator.cpp" -std=c++17 -o testCreator.cpp
g++ "${pathing}/specificTestCreator.cpp" -std=c++17 -o specificTestCreator.cpp


# =========================
# Run pipeline
# =========================

echo "Running Findissues..."
./Findissues.cpp "$SLURM_OUTPUT"
# produces error_info.txt

echo "Running dataCompare..."
./dataCompare.cpp "$ORIGINAL_DATASET" "$MODIFYING_VALUES"
# produces dataset_info.txt and epsilon_and_weights.txt

echo "Running OutputVsData..."
./OutputVsData.cpp "$SLURM_OUTPUT" "$DATASET_INFO"
# this one may need to change, 
#   this code will spew out all specific issues relating
#   to a specific error case. look for targetIndex if 
#   you wish to change the specific error cases we are finding
# produces errorFound.txt

echo "Running testCreator..."
./testCreator.cpp "$SLURM_OUTPUT" "$ORIGINAL_DATASET" "$DATASET_INFO"
# produces UltimateGenericName.hh

echo "Running specificTestCreator..."
./specificTestCreator.cpp "$ULTIMATE_HEADER" "$mode"
# produces MaximumGenericName
#the mode will tell you what do you want to input
#   there is find the largest dsm (you will need to input how many you wish to find)
#   find a specific case
#   find multiple of a specific case

echo "Done."