#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>

using namespace std;

const int TAKE_ROWS_PER_CASE = 6;
const int SKIP_ROWS_PER_CASE = 7;
const int ROW_CYCLE = TAKE_ROWS_PER_CASE + SKIP_ROWS_PER_CASE;
const int DENSITY_COLS = 1650;


/*
epsilon
weights

density objects are the last 5 rows

*/

vector<string> splitCSV(const string& line)
{
    vector<string> result;
    string item;
    stringstream ss(line);

    while (getline(ss, item, ','))
        result.push_back(item);

    return result;
}

bool toDouble(const string& text, double& value)
{
    char* end = nullptr;
    value = strtod(text.c_str(), &end);
    return end != text.c_str();
}

vector<vector<double>> readDensityCases(const string& inputFile)
{
    ifstream file(inputFile);

    if (!file)
    {
        cerr << "Could not open CSV file: " << inputFile << endl;
        return {};
    }

    vector<vector<double>> selectedRows;
    string line;

    // Skip the first/header row
    getline(file, line);

    int dataRowIndex = 0;

    while (getline(file, line))
    {
        if (line.empty())
            continue;

        int cyclePosition = dataRowIndex % 13;

        // Take rows 0-5 of each 13-row block
        // This means CSV rows: 2,3,4,5,6,7 then skip 8-14
        if (cyclePosition < 6)
        {
            vector<string> cells = splitCSV(line);

            if (cells.size() <= 2)
            {
                dataRowIndex++;
                continue;
            }

            vector<double> rowValues;
            bool validRow = true;

            for (size_t i = 2; i < cells.size(); i++)
            {
                double value;

                if (!toDouble(cells[i], value))
                {
                    validRow = false;
                    break;
                }

                rowValues.push_back(value);
            }

            if (validRow)
                selectedRows.push_back(rowValues);
        }

        dataRowIndex++;
    }

    return selectedRows;
}
void writeHeaderFile(
    const string& outputFile,
    const vector<vector<double>>& densityCases
)
{
    ofstream out(outputFile);

    if (!out)
    {
        cerr << "Could not create header file: " << outputFile << endl;
        return;
    }

    int num_cases = static_cast<int>(densityCases.size());

    out << "#ifndef EXTRAPOLATION_HH\n";
    out << "#define EXTRAPOLATION_HH\n\n";

    out << "const int Num_cases = " << num_cases << ";\n\n";
    //out << "const int Density_cols = " << DENSITY_COLS << ";\n\n";

    out << "double density_RK[" << num_cases << "][" << DENSITY_COLS << "] = {\n";

    for (int i = 0; i < num_cases; i++)
    {
        out << "    { ";

        for (int j = 0; j < DENSITY_COLS; j++)
        {
            out << setprecision(17) << densityCases[i][j];

            if (j < DENSITY_COLS - 1)
                out << ", ";
        }

        out << " }";

        if (i < num_cases - 1)
            out << ",";

        out << "\n";
    }

    out << "};\n\n";
    out << "#endif\n";

    cout << "Header file created: " << outputFile << endl;
    cout << "Number of cases found: " << num_cases << endl;
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "just need 1 file, and that is the UltimateGenericName.hh file\n";
        return 1;
    }
    string RK = argv[1];

    string inputFile =
        "./solvingProblems/" + RK;

    string outputDirectory =
    "./GeneratedFiles/";

    string outputFile = outputDirectory + "test_extrapolation.hh";

    vector<vector<double>> densityCases = readDensityCases(inputFile);

    if (densityCases.empty())
    {
        cerr << "No valid cases found." << endl;
        return 1;
    }

    // Write the header file
    writeHeaderFile(outputFile, densityCases);

    cout << "Finished successfully!" << endl;

    return 0;
}
