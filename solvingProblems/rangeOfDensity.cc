#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

/*
    creates a header file from selected rows of a CSV data file.

    required inputs:
        1. CSV data file
        2. Beginning line number
        3. Final line number

    optional arguments:
        --output filename.hh
        --add 2,15,20
        --exclude 6,8

    example:
        rangeDensityCases.sh input.csv 4 10 --add 2,15 --delete 7 --output DensityRange.hh
        rangeDensityCases.sh input.csv 6 15 --add 1,20 --delete 15,13,11,9
        rangeDensityCases.sh input.csv 6 15 --add 1,20,23,25,26
        rangeDensityCases.sh input.csv 6 15
        
        must run it with the shell script file (rangeDensityCases). 
        otherwise the inputs will not work 

    line numbers are 1-based:
        line 1 = first non-empty CSV row

    generated arrays:
        Num_cases
        source_line_cases
        x_cases
        dx_cases
        density_cases
    
        why do i have this? mostly just incase we need it. 
        its unlikely we will need the delete and add function.
        but it is useful incase we do need it.
*/

struct DensityCase
{
    int sourceLineNumber = 0;
    string x;
    string dx;
    vector<string> density;
};

string trim(const string& text)
{
    const size_t first = text.find_first_not_of(" \t\r\n");

    if (first == string::npos)
    {
        return "";
    }

    const size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

vector<string> splitCSV(const string& line)
{
    vector<string> values;
    string current;
    bool insideQuotes = false;

    for (size_t i = 0; i < line.size(); ++i)
    {
        const char c = line[i];

        if (c == '"')
        {
            if (insideQuotes && i + 1 < line.size() && line[i + 1] == '"')
            {
                current += '"';
                ++i;
            }
            else
            {
                insideQuotes = !insideQuotes;
            }
        }
        else if (c == ',' && !insideQuotes)
        {
            values.push_back(trim(current));
            current.clear();
        }
        else
        {
            current += c;
        }
    }

    values.push_back(trim(current));
    return values;
}

vector<int> parseLineList(const string& text)
{
    vector<int> result;
    string item;
    stringstream input(text);

    while (getline(input, item, ','))
    {
        item = trim(item);

        if (item.empty())
        {
            continue;
        }

        const int lineNumber = stoi(item);

        if (lineNumber <= 0)
        {
            throw invalid_argument("Line numbers must be greater than zero.");
        }

        result.push_back(lineNumber);
    }

    return result;
}

vector<string> loadNonEmptyCSVLines(const string& filename)
{
    ifstream input(filename);

    if (!input)
    {
        throw runtime_error("Could not open CSV file: " + filename);
    }

    vector<string> lines;
    string line;

    while (getline(input, line))
    {
        if (!trim(line).empty())
        {
            lines.push_back(line);
        }
    }

    return lines;
}

DensityCase parseDensityCase(const string& csvLine, int sourceLineNumber)
{
    const vector<string> values = splitCSV(csvLine);

    if (values.size() < 3)
    {
        throw runtime_error(
            "CSV line " + to_string(sourceLineNumber) +
            " has fewer than 3 columns. Expected x, dx, and density values."
        );
    }

    DensityCase result;
    result.sourceLineNumber = sourceLineNumber;
    result.x = values[0];
    result.dx = values[1];

    for (size_t column = 2; column < values.size(); ++column)
    {
        result.density.push_back(values[column]);
    }

    return result;
}

void writeHeader(
    const string& outputFilename,
    const vector<DensityCase>& cases
)
{
    if (cases.empty())
    {
        throw runtime_error("No density cases were selected.");
    }

    const size_t densitySize = cases.front().density.size();

    for (const DensityCase& currentCase : cases)
    {
        if (currentCase.density.size() != densitySize)
        {
            throw runtime_error(
                "CSV line " + to_string(currentCase.sourceLineNumber) +
                " has a different number of density values."
            );
        }
    }

    ofstream output(outputFilename);

    if (!output)
    {
        throw runtime_error("Could not create header file: " + outputFilename);
    }

    output << "#ifndef DENSITY_RANGE_CASES_HH\n";
    output << "#define DENSITY_RANGE_CASES_HH\n\n";

    output << "#include <cstddef>\n\n";

    output << "constexpr std::size_t Num_cases = "
           << cases.size() << ";\n";

    output << "constexpr std::size_t Density_size = "
           << densitySize << ";\n\n";

    output << "int source_line_cases[Num_cases] = {\n";

    for (size_t i = 0; i < cases.size(); ++i)
    {
        output << "    " << cases[i].sourceLineNumber;

        if (i + 1 < cases.size())
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n\n";

    output << "double x_cases[Num_cases] = {\n";

    for (size_t i = 0; i < cases.size(); ++i)
    {
        output << "    " << cases[i].x;

        if (i + 1 < cases.size())
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n\n";

    output << "double dx_cases[Num_cases] = {\n";

    for (size_t i = 0; i < cases.size(); ++i)
    {
        output << "    " << cases[i].dx;

        if (i + 1 < cases.size())
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n\n";

    output << "double density_cases[Num_cases][Density_size] = {\n";

    for (size_t i = 0; i < cases.size(); ++i)
    {
        output << "    { ";

        for (size_t j = 0; j < cases[i].density.size(); ++j)
        {
            output << cases[i].density[j];

            if (j + 1 < cases[i].density.size())
            {
                output << ", ";
            }
        }

        output << " }";

        if (i + 1 < cases.size())
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n\n";
    output << "#endif\n";
}

void printUsage(const char* programName)
{
    cerr
        << "Usage:\n"
        << "    " << programName
        << " input.csv beginning_line final_line"
        << " [--add line1,line2,...]"
        << " [--exclude line1,line2,...]"
        << " [--output output.hh]\n\n"
        << "Examples:\n"
        << "    " << programName
        << " data.csv 4 10\n\n"
        << "    " << programName
        << " data.csv 4 10 --add 2,15 --exclude 7\n\n"
        << "    " << programName
        << " data.csv 4 10 --output DensityRange.hh\n";
}

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        printUsage(argv[0]);
        return 1;
    }

    try
    {
        const string inputFilename = argv[1];
        const int beginningLine = stoi(argv[2]);
        const int finalLine = stoi(argv[3]);

        if (beginningLine <= 0 || finalLine <= 0)
        {
            throw invalid_argument("Beginning and final lines must be greater than zero.");
        }

        if (beginningLine > finalLine)
        {
            throw invalid_argument(
                "Beginning line cannot be greater than final line."
            );
        }

        string outputFilename = "DensityRangeCases.hh";
        vector<int> addedLines;
        vector<int> excludedLines;

        for (int argument = 4; argument < argc; ++argument)
        {
            const string option = argv[argument];

            if (option == "--add")
            {
                if (argument + 1 >= argc)
                {
                    throw invalid_argument("--add requires a comma-separated line list.");
                }

                addedLines = parseLineList(argv[++argument]);
            }
            else if (option == "--delete")
            {
                if (argument + 1 >= argc)
                {
                    throw invalid_argument(
                        "--exclude requires a comma-separated line list."
                    );
                }

                excludedLines = parseLineList(argv[++argument]);
            }
            else if (option == "--output")
            {
                if (argument + 1 >= argc)
                {
                    throw invalid_argument("--output requires a filename.");
                }

                outputFilename = argv[++argument];
            }
            else
            {
                throw invalid_argument("Unknown option: " + option);
            }
        }

        const vector<string> csvLines = loadNonEmptyCSVLines(inputFilename);

        if (csvLines.empty())
        {
            throw runtime_error("The CSV file contains no non-empty rows.");
        }

        set<int> selectedLineNumbers;

        for (int lineNumber = beginningLine;
             lineNumber <= finalLine;
             ++lineNumber)
        {
            selectedLineNumbers.insert(lineNumber);
        }

        for (const int lineNumber : addedLines)
        {
            selectedLineNumbers.insert(lineNumber);
        }

        for (const int lineNumber : excludedLines)
        {
            selectedLineNumbers.erase(lineNumber);
        }

        vector<DensityCase> selectedCases;

        for (const int lineNumber : selectedLineNumbers)
        {
            if (lineNumber > static_cast<int>(csvLines.size()))
            {
                throw out_of_range(
                    "Requested line " + to_string(lineNumber) +
                    ", but the CSV only has " +
                    to_string(csvLines.size()) + " non-empty rows."
                );
            }

            selectedCases.push_back(
                parseDensityCase(csvLines[lineNumber - 1], lineNumber)
            );
        }

        writeHeader(outputFilename, selectedCases);

        cout << "Input CSV: " << inputFilename << "\n";
        cout << "Beginning line: " << beginningLine << "\n";
        cout << "Final line: " << finalLine << "\n";
        cout << "Cases written: " << selectedCases.size() << "\n";
        cout << "Header created: " << outputFilename << "\n\n";

        cout << "Selected CSV lines:\n";

        for (const DensityCase& currentCase : selectedCases)
        {
            cout << "    " << currentCase.sourceLineNumber << "\n";
        }

        return 0;
    }
    catch (const exception& error)
    {
        cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
/*
int main2(int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr
            << "Usage:\n"
            << argv[0]
            << " input.csv beginning_line ending_line"
            << " [--add lines] [--delete lines]\n";

        return 1;
    }

    string inputFile = argv[1];
    int beginningLine = stoi(argv[2]);
    int endingLine = stoi(argv[3]);

    string addCases = "";
    string deleteCases = "";

    for (int i = 4; i < argc; i++)
    {
        string option = argv[i];

        if (option == "--add")
        {
            if (i + 1 >= argc)
            {
                cerr << "--add must be followed by line numbers.\n";
                return 1;
            }

            addCases = argv[++i];
        }
        else if (option == "--delete" || option == "--exclude")
        {
            if (i + 1 >= argc)
            {
                cerr << option
                     << " must be followed by line numbers.\n";

                return 1;
            }

            deleteCases = argv[++i];
        }
        else
        {
            cerr << "Unknown option: " << option << "\n";
            return 1;
        }
    }

    cout << "Input file: " << inputFile << "\n";
    cout << "Beginning line: " << beginningLine << "\n";
    cout << "Ending line: " << endingLine << "\n";

    if (!addCases.empty())
    {
        cout << "Cases to add: " << addCases << "\n";
    }

    if (!deleteCases.empty())
    {
        cout << "Cases to delete: " << deleteCases << "\n";
    }
    return 0;
}
    */

/*
find a range of cases specified in a csv file. say row / density cases 4 to 10.
    should print out all density cases between 4 through 10. 
an added feature we should have, just incase, is adding specific cases.
    both adding multiple specific cases, and/or excluding a specific case
    within a range of code. 


*/
