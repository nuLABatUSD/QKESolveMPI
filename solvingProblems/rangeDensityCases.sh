#C:/msys64/usr/bin/bash.exe
#!/bin/sh


rangepath="./solvingProblems/rangeOfDensity.cc"
rangeexepath="./solvingProblems/rangeOfDensity.exe"

inputfile="$1"
beginingValue="$2"
finalValue="$3"

g++ "$rangepath" -std=c++17 -O2 -o "$rangeexepath"

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "creating file .... DensityRangeCases.hh"
: '
"$rangeexepath" \
    "$inputfile" \
    "$beginingValue" \
    "$finalValue"
 '

"$rangeexepath" "$@"

: '
#C:/msys64/usr/bin/bash.exe
#!/bin/bash

compiler="C:/msys64/mingw64/bin/g++.exe"

source="./solvingProblems/rangeOfDensity.cpp"
program="./solvingProblems/rangeOfDensity.exe"

echo "Compiling..."

"$compiler" "$source" -std=c++17 -O2 -o "$program"

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "Running..."

"$program" "$@"
'
#idk why there are multiple different types of commend such as ''' ''' and : ' ' 
#   and << 'comment' ... comment
#   but what ever works 