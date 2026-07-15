#C:/msys64/usr/bin/bash.exe
#!/bin/bash


pathing="./solvingProblems"

#lastLine="59-QKE_fix_mat.csv"
#RKoutput (2).csv
read -p "enter lastLine file name: " lastLine



g++ ./solvingProblems/restartCreator.cpp -std=c++17 -o restartCreator.cpp

echo "running restartCreator.cpp"

./restartCreator.cpp "$pathing/$lastLine"