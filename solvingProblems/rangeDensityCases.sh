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

#idk why there are multiple different types of commend such as ''' ''' and : ' ' 
#   and << 'comment' ... comment
#   but what ever works 
