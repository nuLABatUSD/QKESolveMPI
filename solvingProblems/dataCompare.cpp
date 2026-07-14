
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <map>
#include <vector>
#include <algorithm>

#include <cmath>
#include <iomanip>
#include <sstream>

//#include "FINDissues.cpp"

using namespace std;

/*
206 bins

- find bins
- tell if its px, py, pz
- match the x in the output to the fix_mat value
    is going to match the first colum of the data file 

e_coll_eps.csv

first row are the epsilon values 
2nd weights 

59-QKE_fix_mat.csv

when we see a lot of steps, we need to be able to find that location with many steps.
    for example that one x value with 2.72735e+15 or another x value with 8.48641e+14
    they have been repeated several times. and find out why they were being repeated 

    - find the most amounts of repeated processes (most reoccuring values for x)
    - find that index number asssociated with it 
    - 

    get output file, SLURM and compare the x value 
    to the data file 59-QKE-fix-mat
*/

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "must have 2 files, dataCompare dataset.csv and epsilon.csv\n";
        return 1;
    }
    string dataset = argv[1];
    string eps = argv[2];

    string file1 = "./solvingProblems/" + eps;
    //location of the epsilon and weight values are

    string file2 = "./solvingProblems/" + dataset;
    //location of the datasheet is 

    regex indexPattern(R"(index\s*=\s*([0-9]+))");

    ifstream f1(file1);
    ifstream f2(file2);

    string line;
    string line2;
    string currentX = "unknown";

    //when dealing with the data and plotting it
    /*
    comparing the 59-qke-fix. each data point is actually a measurement in 4
    o, x, y, z; 
    basically every 1st positino is x, with it ending with dx, y4, then y5 
    so its 
    x, dx, y4, y5
    from here we can compare that to the weighted values from the other data sheet
    06-nu-e-coll whatever. 

    */
    //ifstream file(filename);
        if (!f2)
        {
            cerr << "Could not open file: " << file2 << endl;
            return 1;
        }
        if (!f1)
        {
            cerr << "Could not open file: " << file2 << endl;
            return 1;
        }


        ofstream out("dataset_info.txt");
        //allows is to create a condensed data file of the data sheet
        //  only looking at the first 2 columns, and from there use that info to compare

        ofstream out2("epsilon_and_weights.txt");
        //allows us to see all of the epsilon and weight values


    if (!out)
    {
        cerr << "Could not create dataset_info.txt\n";
        return 1;
    }
    if (!out2)
    {
        cerr << "Could not create epsilon_and_weights.txt\n";
        return 1;
    }

    vector<string> xValues;
    vector<string> dxValues;
    string cell;
    long long totalValues = 0;

    while (getline(f2, line))
    {
        stringstream ss(line);

        size_t colIndex = 0;

        while (getline(ss, cell, ','))
        {
            if (cell.empty())
            continue;

        if (colIndex == 0)
        {
            xValues.push_back(cell);
        }
        else if (colIndex == 1)
        {
            dxValues.push_back(cell);
            break; // done with this row, move to next row
        }

        colIndex++;
    }


    }


    //long long totalDatasets = totalValues / 4;

    out << "Total values = " << totalValues << "\n";

    out << "X values:\n";

    for (size_t i = 0; i < xValues.size(); i++)
    {
        out << std::left
         << "x["<< std::setw(3) << i <<"] = "
         << std::setw(22) << xValues[i]   // make room for long doubles

         << "\tdx["<< std::setw(3) << i <<"] = "
         << std::setw(22) << dxValues[i]

         << '\n' << '\n';

    }

    if (totalValues % 4 != 0)
    {
        out << "\nWarning: leftover values = " << (totalValues % 4) << "\n";
    }

    f2.close();
    out.close();
    cout << "Done. Results written to dataset_info.txt.txt\n";

    vector<string> EWValues;
    string cell2;
    while (getline(f1, line2))
    {
        stringstream ss2(line2);

        while (getline(ss2, cell2, ','))
        {
            if (cell2.empty()){
                continue;
            }
            else{
                EWValues.push_back(cell2);
            }
        }
    }

    for (size_t i = 0; i < 206; i++)
    {
        out2 << std::left
         << "epsilon[" << std::setw(3) << i <<"] = "
         << std::setw(22) << EWValues[i]   // make room for long doubles
         << "weights[" << (i) << "] = "
         << EWValues[206+i]
         << '\n';


    }

    cout << "Done. Results written to epsilon_and_weights.txt\n";

    f1.close();
    out2.close();

    return 0;
}

struct FileShape
{
    size_t rows = 0;
    size_t cols = 0;
};

FileShape GetFileShape(const std::string& filename)
{
    std::ifstream file(filename);
    FileShape shape;

    if (!file)
    {
        std::cerr << "Failed to open " << filename << "\n";
        return shape;
    }

    std::string line;

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        shape.rows++;

        std::stringstream ss(line);
        std::string cell;
        size_t currentCols = 0;

        while (std::getline(ss, cell, ','))
        {
            if (cell.empty())
                continue;

            currentCols++;
        }

        if (currentCols > shape.cols)
            shape.cols = currentCols;
    }

    return shape;
}