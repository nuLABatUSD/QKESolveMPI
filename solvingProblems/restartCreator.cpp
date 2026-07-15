#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <filesystem>

using namespace std;

const int DENSITY_SIZE = 1650;
const int EXPECTED_VALUES = DENSITY_SIZE + 2;


/*
Splits one CSV row at every comma.
*/
vector<string> splitCSV(const string& line)
{
    vector<string> result;
    string item;
    stringstream ss(line);

    while (getline(ss, item, ','))
    {
        result.push_back(item);
    }

    return result;
}


/*
Reads the CSV and returns only its last nonempty row.
*/
string getLastCSVRow(const string& csvFileName)
{
    ifstream file(csvFileName);

    if (!file.is_open()) {
        std::cout << "Error: Could not open file!"<<endl;;
    }


    if (!file)
    {
        cout << "Trying to open: "
         << filesystem::absolute(csvFileName) << endl;
        cerr << "Could not open CSV file: "
             << csvFileName << endl;

        return "";
    }

    string line;
    string lastRow;

    while (getline(file, line))
    {
        // Check whether the line contains any actual characters
        // besides spaces, tabs, or Windows carriage returns.
        if (line.find_first_not_of(" \t\r\n") != string::npos)
        {
            lastRow = line;
        }
    }

    file.close();

    return lastRow;
}


/*
Creates the restart header file.
*/
bool writeRestartHeader(
    const string& headerFileName,
    const vector<string>& values
)
{
    
    
    if (values.size() < EXPECTED_VALUES)
    {
        cerr << "The final CSV row only has "
             << values.size()
             << " values." << endl;

        cerr << "Expected at least "
             << EXPECTED_VALUES
             << " values: x, dx, and "
             << DENSITY_SIZE
             << " density values." << endl;

        return false;
    }
        

    ofstream out(headerFileName);

    if (!out)
    {
        cerr << "Could not create header file: "
             << headerFileName << endl;

        return false;
    }

    out << "#ifndef RESTART_VALUES_HH\n";
    out << "#define RESTART_VALUES_HH\n\n";

    /*
    The strings from the CSV are written directly into
    the header so their original precision is preserved.
    */
    out << "double x_restart = "
        << values[0]
        << ";\n\n";

    out << "double dx_restart = "
        << values[1]
        << ";\n\n";

    out << "double density_restart["
        << DENSITY_SIZE
        << "] = {\n";

    for (int i = 0; i < DENSITY_SIZE; i++)
    {
        out << "    " << values[i + 2];

        if (i < DENSITY_SIZE - 1)
        {
            out << ",";
        }

        out << "\n";
    }

    out << "};\n\n";
    out << "#endif\n";

    out.close();

    cout << "Header file created: "
         << headerFileName << endl;

    return true;
}


int main(int argc, char* argv[])
{
    
    if (argc != 2)
    {
        cerr << "Usage: restartCreator <CSV filename>\n";
        return 1;
    }

    //59-QKE_fix_mat
    //RKoutput (2).csv
    string csvName = argv[1];

    string csvFile =
        "./" + csvName;

    string headerFile =
        "./solvingProblems/restart.hh";

    cout << "Current working directory: "
         << filesystem::current_path() << endl;

    cout << "Trying to open: "
         << filesystem::absolute(csvFile) << endl;

    string lastRow = getLastCSVRow(csvFile);

    if (lastRow.empty())
    {
        cerr << "The CSV file is empty or could not be read.\n";
        return 1;
    }

    vector<string> values = splitCSV(lastRow);

    if (!writeRestartHeader(headerFile, values))
    {
        return 1;
    }

    return 0;
}
