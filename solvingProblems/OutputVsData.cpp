
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

#include <chrono>


using namespace std;

struct ProblemCase
{
    double dsm;
    double x;
    double dx;
    std::string valueLine;
    std::string problemLine;
    int fileALineNumber;
};
struct CondensedRow
{
    int lineNumber;
    int index;
    double x;
    double dx;
    std::string fullLine;
};

string sciKey(double value, int precision = 4)
{
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(precision) << value;
    return ss.str();
}

chrono::high_resolution_clock::time_point StartTimer()
{
    return std::chrono::high_resolution_clock::now();
}

std::string xSearchKey(double value)
{
    if (value == 0.0)
        return "0";

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(15) << value;

    std::string s = ss.str();

    size_t ePos = s.find('e');
    if (ePos != std::string::npos)
        s = s.substr(0, ePos);

    s.erase(std::remove(s.begin(), s.end(), '.'), s.end());
    s.erase(std::remove(s.begin(), s.end(), '-'), s.end());

    while (!s.empty() && s.back() == '0')
        s.pop_back();

    if (s.size() > 1)
        s.pop_back();

    return s;
}

void EndTimer(const std::chrono::high_resolution_clock::time_point& start,
              const std::string& label = "Elapsed time")
{
    auto end = std::chrono::high_resolution_clock::now();

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::chrono::duration<double> sec = end - start;

    std::cout << label << ": "
              << ms.count() << " ms ("
              << sec.count() << " s)\n";
}

bool extractValues(const std::string& line, double& dsm, double& x, double& dx)
{
    std::regex valuePattern(
        R"(dsm\s*=\s*([-+0-9.eE]+),\s*x\s*=\s*([-+0-9.eE]+),\s*dx\s*=\s*([-+0-9.eE]+))"
    );

    std::smatch match;

    if (std::regex_search(line, match, valuePattern))
    {
        dsm = std::stod(match[1]);
        x   = std::stod(match[2]);
        dx  = std::stod(match[3]);
        return true;
    }

    return false;
}

std::string mantissaKey(double value, int decimals = 4)
{
    if (value == 0.0)
        return "0.0000";

    double exponent = std::floor(std::log10(std::abs(value)));
    double mantissa = value / std::pow(10.0, exponent);

    double factor = std::pow(10.0, decimals);

    // truncate, do not round
    mantissa = std::trunc(mantissa * factor) / factor;

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(decimals) << mantissa;

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

string roundedSciKey(double value, int precision = 5)
{
    ostringstream ss;
    ss << scientific << setprecision(precision) << value;
    return ss.str();
}

int findMatchingLineInFileB(
    const std::string& fileBName,
    double targetX,
    std::string& matchedLine
)
{
    std::ifstream fileB(fileBName);
    if (!fileB) return -1;

    std::string targetXKey = xSearchKey(targetX);

    std::regex xPattern(R"(x\s*=\s*([-+0-9.eE]+))");

    std::string line;
    int lineNumber = 0;

    while (std::getline(fileB, line))
    {
        lineNumber++;

        std::smatch xMatch;

        if (std::regex_search(line, xMatch, xPattern))
        {
            double fileX = std::stod(xMatch[1]);
            std::string fileXKey = xSearchKey(fileX);

            if (fileXKey == targetXKey)
            {
                matchedLine = line;
                return lineNumber;
            }
        }
    }

    return -1;
}

std::string expandScientific(double value)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(0) << value;
    return ss.str();
}

std::string convertExpandedX(const std::string& expanded)
{
    std::string result;

    for (char c : expanded)
    {
        if (std::isdigit(c))
            result += c;
    }

    // drop last digit to avoid rounding issue
    if (result.size() > 1)
        result.pop_back();

    return result;
}

void DebugSearchForX(const vector<CondensedRow>& rows, double targetX, int maxPrint = 20)
{
    string targetKey = roundedSciKey(targetX, 5);

    cout << "\n================ DEBUG SEARCH FOR X ================\n";
    cout << "Target original x = " << setprecision(17) << targetX << "\n";
    cout << "Target sci key    = " << targetKey << "\n\n";

    int printed = 0;

    for (const CondensedRow& row : rows)
    {
        string rowKey = roundedSciKey(row.x, 5);

        bool match = rowKey == targetKey;

        if (match || printed < maxPrint)
        {
            cout << "Row line " << row.lineNumber
                 << " | index " << row.index << "\n";
            cout << "  row original x = " << setprecision(17) << row.x << "\n";
            cout << "  row sci key    = " << rowKey << "\n";
            cout << "  compare        = " << targetKey << " vs " << rowKey << "\n";
            cout << "  match?         = " << (match ? "YES" : "NO") << "\n";
            cout << "  full line      = " << row.fullLine << "\n";
            cout << "----------------------------------------------------\n";

            printed++;
        }

        if (match)
            break;
    }

    cout << "====================================================\n";
}

bool roughlyEqual(double a, double b, double relativeTolerance = 1e-4)
{
    double diff = std::abs(a - b);
    double scale = std::max(std::abs(a), std::abs(b));

    if (scale == 0.0)
        return diff < relativeTolerance;

    return diff <= relativeTolerance * scale;
}

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

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "must have 2 inputs, a slurm file and the condensed data set (dataset_info)\n";
        return 1;
    }
    string slurm = argv[1];
    string dataset_info = argv[2];

    //slurm-18329938-Copy   output data 1
    //slurm-18475761        output data 2
    string fileAName = "./solvingProblems/" + slurm;
    //string fileBName = "C:/Users/astar/Desktop/workspace/QKESolveMPI-main/solvingProblems/59-QKE_fix_mat.csv";
    string fileBName = "./solvingProblems/" + dataset_info;

    ifstream fileA(fileAName);
    ofstream out("errorFound.txt");

    std::vector<CondensedRow> condensedRows = loadCondensedFile(fileBName);
    cout << "Condensed rows loaded: " << condensedRows.size() << endl;

    if (condensedRows.empty())
    {
        cerr << "ERROR: condensedRows is empty. Regex probably did not match the file format, or file path is wrong.\n";
        return 1;
    }

    auto start = std::chrono::high_resolution_clock::now();

    if (!fileA)
    {
        std::cerr << "Failed to open " << fileAName << "\n";
        return 1;
    }

    if (!out)
    {
        std::cerr << "Failed to open errorFound.txt\n";
        return 1;
    }

    const int targetIndex = 818;
    /*
    index values = 818 appeared 178 times		2 is the index value thus P_y is the problem
    index values = 814 appeared 174 times		2 is the index value thus P_y is the problem
    index values = 810 appeared 166 times		2 is the index value thus P_y is the problem
    index values = 806 appeared 162 times		2 is the index value thus P_y is the problem
    index values = 822 appeared 138 times		2 is the index value thus P_y is
    */
   /*
   818 
   814 
   810 
   806 
   822 
   821 has 1 data point but the largest dsm value


   */

    std::regex indexPattern(R"(problem index\s*=\s*(\d+))");

    std::string line;
    std::string previousLine;
    int lineNumber = 0;
    int casesFound = 0;

    

    while (std::getline(fileA, line))
    {
        lineNumber++;

        std::smatch indexMatch;

        if (std::regex_search(line, indexMatch, indexPattern))
        {
            int foundIndex = std::stoi(indexMatch[1]);

            if (foundIndex == targetIndex)
            {
                double dsm, x, dx;

                if (extractValues(previousLine, dsm, x, dx))
                {
                    casesFound++;

                    CondensedRow matchedRow;
                    bool foundMatch = findMatchingCondensedRow(condensedRows, x, matchedRow);

                    //int fileBLine = findMatchingLineInFileB(fileBName, x, matchedLine);

                    out << "========================================\n";
                    out << "Problem case #" << casesFound << "\n";
                    out << "Found in | slurm-18329938-Copy.txt | at line: " << lineNumber << "\n";
                    out << previousLine << "\n";
                    out << line << "\n\n";

                    //out << "Extracted values:\n";
                    //out << "dsm = " << dsm << "\n";
                    //out << "x   = " << x << "\n";
                    //out << "dx  = " << dx << "\n\n";

                    if (!condensedRows.empty())
                    {
                        DebugSearchForX(condensedRows, x, 10);
                    }
                    //ignore this, this is a debugging statement

                    if (foundMatch)
                    {
                        out << "Matching x found in condensed file on line: "
                            << matchedRow.lineNumber << "\n";

                        out << "Condensed index: " << matchedRow.index << "\n";
                        out << matchedRow.fullLine << "\n";
                    }
                    else
                    {
                        out << "No matching x found in condensed file.\n";
                    }

                    out << "\n";
                }
                else
                {
                    out << "========================================\n";
                    out << "Problem index found at | slurm | line "
                        << lineNumber << ", but could not read dsm/x/dx from the line above.\n";
                    out << "Line above was:\n";
                    out << previousLine << "\n\n";
                }
            }
        }

        previousLine = line;
    }

    out << "========================================\n";
    out << "Total problem index cases found: " << casesFound << "\n";

    fileA.close();
    out.close();
    

    cout << "=======================================================" << endl;
    cout << "the task is finished" << endl;
    cout << "=======================================================" << endl;

    EndTimer(start, "Loop runtime"); // end timer + print result
    
    return 0;
}
