#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <regex>
#include <unordered_map>
#include <iomanip>
#include <unordered_set>


#include <cmath>
#include <algorithm>

using namespace std;
/*
so dont create a file that will plot all possible error cases
    grab 1 error case
    grab its x value
    grab the first dx value (if there are mutliple x error cases)
    and its density variables
    put all of those issues on 1 test case, it should have 1 element
        per array.
    
*/

struct CaseData
{
    double x;
    double dx;
    double dsm;
    int problemIndex;
    vector<double> density;
};
struct CondensedRow
{
    int lineNumber;
    int index;
    double x;
    double dx;
    std::string fullLine;
};

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

vector<string> loadDataSheetLines(const string& dataFile)
{
    vector<string> lines;

    ifstream file(dataFile);

    if (!file)
    {
        cerr << "Could not open data sheet: " << dataFile << endl;
        return lines;
    }

    string line;

    while (getline(file, line))
    {
        if (!line.empty())
            lines.push_back(line);
    }

    return lines;
}
string roundedSciKey(double value, int precision = 5)
{
    ostringstream ss;
    ss << scientific << setprecision(precision) << value;
    return ss.str();
}
std::vector<CondensedRow> loadCondensedFile(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<CondensedRow> rows;

    if (!file)
    {
        std::cerr << "Failed to open " << filename << "\n";
        return rows;
    }

    std::regex pattern(
        R"(x\[\s*(\d+)\s*\]\s*=\s*([-+0-9.eE]+)\s*dx\[\s*\d+\s*\]\s*=\s*([-+0-9.eE]+))"
    );

    std::string line;
    int lineNumber = 0;

    while (std::getline(file, line))
    {
        lineNumber++;

        if (line.empty())
            continue;

        std::smatch match;

        if (std::regex_search(line, match, pattern))
        {
            CondensedRow row;
            row.lineNumber = lineNumber;
            row.index = std::stoi(match[1]);
            row.x = std::stod(match[2]);
            row.dx = std::stod(match[3]);
            row.fullLine = line;

            rows.push_back(row);
        }
    }

    return rows;
}
bool extractProblemIndex(const string& line, int& problemIndex)
{
    regex pattern(R"(problem index\s*=\s*(\d+))");
    smatch match;

    if (regex_search(line, match, pattern))
    {
        problemIndex = stoi(match[1]);
        return true;
    }

    return false;
}
int findDataRowByX(const vector<string>& dataLines, double targetX)
{
    string targetKey = roundedSciKey(targetX, 5);

    for (int i = 0; i < static_cast<int>(dataLines.size()); i++)
    {
        vector<string> values = splitCSV(dataLines[i]);

        if (values.empty())
            continue;

        double rowX = stod(values[0]);
        string rowKey = roundedSciKey(rowX, 5);

        if (rowKey == targetKey)
            return i;
    }

    return -1;
}
bool extractProblemValues(const string& line, double& dsm, double& x, double& dx)
{
    regex pattern(R"(dsm\s*=\s*([-+0-9.eE]+)\s*,\s*x\s*=\s*([-+0-9.eE]+)\s*,\s*dx\s*=\s*([-+0-9.eE]+))");

    smatch match;

    if (regex_search(line, match, pattern))
    {
        dsm = stod(match[1]);
        x   = stod(match[2]);
        dx  = stod(match[3]);
        return true;
    }
    return false;
}
//a helper function

bool findMatchingCondensedRow(
    const vector<CondensedRow>& rows,
    double targetX,
    CondensedRow& matchedRow
)
{
    string targetKey = roundedSciKey(targetX, 5);

    for (const CondensedRow& row : rows)
    {
        string rowKey = roundedSciKey(row.x, 5);

        if (rowKey == targetKey)
        {
            matchedRow = row;
            return true;
        }
    }

    return false;
}
void processProblemFile2(
    const string& problemFile,
    const vector<string>& dataLines,
    vector<CaseData>& cases
)
{
    ifstream file(problemFile);

    if (!file)
    {
        cerr << "Could not open problem file: " << problemFile << endl;
        return;
    }

    regex condensedIndexPattern(R"(Condensed index:\s*(\d+))");

    unordered_set<int> usedCondensedIndexes;

    string line;

    double currentProblemDx = 0.0;
    double currentProblemDsm = 0.0;
    double currentProblemX = 0.0;

    bool hasProblemValues = false;

    while (getline(file, line))
    {
        double tempDsm;
        double tempX;
        double tempDx;

        if (extractProblemValues(line, tempDsm, tempX, tempDx))
        {
            currentProblemDsm = tempDsm;
            currentProblemX = tempX;
            currentProblemDx = tempDx;
            hasProblemValues = true;

            cout << "Found problem values: "
                 << "dsm = " << currentProblemDsm
                 << ", x = " << currentProblemX
                 << ", dx = " << currentProblemDx << endl;
        }

        smatch match;

        if (regex_search(line, match, condensedIndexPattern))
        {
            int condensedIndex = stoi(match[1]);

            if (usedCondensedIndexes.count(condensedIndex) > 0)
            {
                cout << "Skipping repeated x / condensed index: "
                     << condensedIndex << endl;

                hasProblemValues = false;
                continue;
            }

            usedCondensedIndexes.insert(condensedIndex);

            int dataRowIndex = condensedIndex;

            if (dataRowIndex < 0 || dataRowIndex >= static_cast<int>(dataLines.size()))
            {
                cout << "Data row index out of range: " << dataRowIndex << endl;
                hasProblemValues = false;
                continue;
            }

            vector<string> values = splitCSV(dataLines[dataRowIndex]);

            if (values.size() < 3)
            {
                cout << "Bad data row at index: " << dataRowIndex << endl;
                hasProblemValues = false;
                continue;
            }

            CaseData newCase;

            newCase.x = stod(values[0]);

            if (hasProblemValues)
            {
                newCase.dx = currentProblemDx;
                newCase.dsm = currentProblemDsm;
            }
            else
            {
                newCase.dx = stod(values[1]);
                newCase.dsm = 0.0;

                cout << "WARNING: condensed index found without dsm line before it.\n";
            }

            for (size_t i = 2; i < values.size(); i++)
            {
                newCase.density.push_back(stod(values[i]));
            }

            cases.push_back(newCase);

            cout << "Stored case:"
                 << " x = " << newCase.x
                 << ", dx = " << newCase.dx
                 << ", dsm = " << newCase.dsm
                 << endl;

            hasProblemValues = false;
        }
    }
}
//old code, may need to delete later

void processProblemFile(
    const string& problemFile,
    const vector<string>& dataLines,
    const vector<CondensedRow>& condensedRows,
    vector<CaseData>& cases
)
{
    ifstream file(problemFile);

    if (!file)
    {
        cerr << "Could not open problem file: " << problemFile << endl;
        return;
    }

    unordered_set<int> usedCondensedIndexes;

    string line;

    double savedDsm = 0.0;
    double savedX = 0.0;
    double savedDx = 0.0;
    bool hasSavedValues = false;

    while (getline(file, line))
    {
        double dsm;
        double x;
        double dx;

        if (extractProblemValues(line, dsm, x, dx))
        {
            savedDsm = dsm;
            savedX = x;
            savedDx = dx;
            hasSavedValues = true;
            continue;
        }

        int problemIndex;

        if (!extractProblemIndex(line, problemIndex))
            continue;

        if (!hasSavedValues)
            continue;

        CondensedRow matchedRow;

        if (!findMatchingCondensedRow(condensedRows, savedX, matchedRow))
        {
            cout << "No condensed match for x = "
                 << setprecision(17) << savedX << endl;

            hasSavedValues = false;
            continue;
        }

        if (usedCondensedIndexes.count(matchedRow.index) > 0)
        {
            hasSavedValues = false;
            continue;
        }

        usedCondensedIndexes.insert(matchedRow.index);

        int dataRowIndex = matchedRow.index;

        if (dataRowIndex < 0 || dataRowIndex >= static_cast<int>(dataLines.size()))
        {
            cout << "Data row index out of range: " << dataRowIndex << endl;
            hasSavedValues = false;
            continue;
        }

        vector<string> values = splitCSV(dataLines[dataRowIndex]);

        CaseData newCase;

        newCase.x = stod(values[0]);
        newCase.dx = savedDx;
        newCase.dsm = savedDsm;
        newCase.problemIndex = problemIndex;

        for (size_t i = 2; i < values.size(); i++)
        {
            newCase.density.push_back(stod(values[i]));
        }

        cases.push_back(newCase);

        cout << "Stored case:"
             << " x = " << setprecision(17) << newCase.x
             << ", dx = " << newCase.dx
             << ", dsm = " << newCase.dsm
             << ", problem index = " << newCase.problemIndex
             << endl;

        hasSavedValues = false;
    }
}

void writeHeaderFile(
    const string& headerFileName,
    const vector<double>& x_cases,
    const vector<double>& dx_cases,
    const vector<double>& dsm_cases,
    const vector<int>& problem_index_cases,
    const vector<vector<double>>& density_cases
)
{
    ofstream out(headerFileName);

    if (!out)
    {
        cerr << "Could not create header file: " << headerFileName << endl;
        return;
    }

    size_t Num_cases = x_cases.size();

    size_t densitySize = 0;

    if (!density_cases.empty())
    {
        densitySize = density_cases[0].size();
    }

    out << "#ifndef ULTIMATE_GENERIC_NAME_HH\n";
    out << "#define ULTIMATE_GENERIC_NAME_HH\n\n";

    out << "#include <cstddef>\n\n";

    //out << "double* test_density;\n\n";

    out << "int Num_cases = " << Num_cases << ";\n\n";

    out << "/*\n";
    out << "problem index corresponding to the x value\n";

    for (size_t i = 0; i < problem_index_cases.size(); i++)
    {
        out << problem_index_cases[i] << "\n";
    }

    out << "*/\n\n";

    out << "double x_cases[" << Num_cases << "] = {\n";

    for (size_t i = 0; i < Num_cases; i++)
    {
        out << "    " << fixed << setprecision(0) << x_cases[i];

        if (i < Num_cases - 1)
            out << ",";

        out << "\n";
    }

    out << "};\n\n";

    out << defaultfloat;
    out << "double dsm_cases[" << Num_cases << "] = {\n";

    for (size_t i = 0; i < Num_cases; i++)
    {
        out << "    " << setprecision(17) << dsm_cases[i];

        if (i < Num_cases - 1)
            out << ",";

        out << "\n";
    }

    out << "};\n\n";

    out << "double dx_cases[" << Num_cases << "] = {\n";

    for (size_t i = 0; i < Num_cases; i++)
    {
        out << "    " << setprecision(17) << dx_cases[i];

        if (i < Num_cases - 1)
            out << ",";

        out << "\n";
    }

    out << "};\n\n";

    out << "double density_cases[" << Num_cases << "][" << densitySize << "] = {\n";

    for (size_t i = 0; i < Num_cases; i++)
    {
        out << "    { ";

        for (size_t j = 0; j < density_cases[i].size(); j++)
        {
            out << setprecision(17) << density_cases[i][j];

            if (j < density_cases[i].size() - 1)
                out << ", ";
        }

        out << " }";

        if (i < Num_cases - 1)
            out << ",";

        out << "\n";
    }

    out << "};\n\n";

    out << "#endif\n";

    out.close();

    cout << "Header file created: " << headerFileName << endl;
}
/*
this is where it creates the header file for the super computer to run
*/

int testCreation(string outputfile, string dataname, string condensedDataName)
{
    //slurm-18475761.txt
    //65-fix_10_small.csv
    //dataset_info.txt
    vector<string> problemFiles =
    {
        "./solvingProblems/"+outputfile
    };
    //locations for the error case files
    //should be able to go through all the files just fine

    string dataSheetFile = "./solvingProblems/" + dataname;
    //data sheet file path location

    string condensedData = "./solvingProblems/" + condensedDataName;
    //condensed data sheet
    
    vector<CaseData> cases;

    vector<string> dataLines = loadDataSheetLines(dataSheetFile);

    vector<CondensedRow> condensedRows = loadCondensedFile(condensedData);

    if (condensedRows.empty())
    {
        cerr << "No condensed rows loaded.\n";
        return 1;
    }
    for (const string& file : problemFiles)
    {
        processProblemFile(file, dataLines, condensedRows, cases);
    }

    vector<double> x_cases;
    vector<double> dx_cases;
    vector<double> dsm_cases;
    vector<int> problem_index_cases;
    vector<vector<double>> density_cases;

    for (const CaseData& c : cases)
    {
        dsm_cases.push_back(c.dsm);

        x_cases.push_back(c.x);
        dx_cases.push_back(c.dx);
        problem_index_cases.push_back(c.problemIndex);
        density_cases.push_back(c.density);
    }

    int num_cases = static_cast<int>(cases.size());

    cout << "Number of cases found: " << num_cases << endl;

    writeHeaderFile(
        "UltimateGenericName.hh",
        x_cases,
        dx_cases,
        dsm_cases,
        problem_index_cases,
        density_cases
        /*
        create dsm_cases array. 
        that contains the first dsm value
        */
    );
    return 0;
}

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        cerr << "must have 3 input files. the slurm output file, the original data csv file, and the condensed data set(dataset_info)\n";
        return 1;
    }
    string slurm = argv[1];  
    string originalSet = argv[2];  
    string dataset_info = argv[3];  

    vector<string> problemFiles =
    {
        
        "./solvingProblems/" + slurm
        //wanted to make this into an array of locations that you can find the issues from.
        //technically you can, but doesnt want to work with shell script. 
        //so reduced it to 1 file

    };
    //locations for the error case files
    //should be able to go through all the files just fine

    string dataSheetFile = "./solvingProblems/" + originalSet;
    //data sheet file path location

    string condensedData = "./solvingProblems/" + dataset_info;
    //condensed data sheet
    
    vector<CaseData> cases;

    vector<string> dataLines = loadDataSheetLines(dataSheetFile);

    vector<CondensedRow> condensedRows = loadCondensedFile(condensedData);

    if (condensedRows.empty())
    {
        cerr << "No condensed rows loaded.\n";
        return 1;
    }
    for (const string& file : problemFiles)
    {
        processProblemFile(file, dataLines, condensedRows, cases);
    }

    vector<double> x_cases;
    vector<double> dx_cases;
    vector<double> dsm_cases;
    vector<int> problem_index_cases;
    vector<vector<double>> density_cases;

    for (const CaseData& c : cases)
    {
        dsm_cases.push_back(c.dsm);

        x_cases.push_back(c.x);
        dx_cases.push_back(c.dx);
        problem_index_cases.push_back(c.problemIndex);
        density_cases.push_back(c.density);
    }

    int num_cases = static_cast<int>(cases.size());

    cout << "Number of cases found: " << num_cases << endl;

    writeHeaderFile(
        "UltimateGenericName.hh",
        x_cases,
        dx_cases,
        dsm_cases,
        problem_index_cases,
        density_cases
        /*
        create dsm_cases array. 
        that contains the first dsm value
        */
    );

    return 0;
}



