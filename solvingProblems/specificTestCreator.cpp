#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>

using namespace std;

const int DEFAULT_LARGEST_COUNT = 5;
string indexValues [4] = {"P_0", "P_x", "P_y", "P_z"};

string readWholeFile(const string& fileName)
{
    ifstream file(fileName);

    if (!file)
    {
        cerr << "Could not open file: " << fileName << endl;
        return "";
    }

    stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

string trim(const string& s)
{
    size_t start = s.find_first_not_of(" \n\r\t");
    size_t end = s.find_last_not_of(" \n\r\t");

    if (start == string::npos)
        return "";

    return s.substr(start, end - start + 1);
}

vector<string> splitTopLevelCommas(const string& text)
{
    vector<string> result;
    string current;
    int braceDepth = 0;

    for (char c : text)
    {
        if (c == '{')
            braceDepth++;
        else if (c == '}')
            braceDepth--;

        if (c == ',' && braceDepth == 0)
        {
            result.push_back(trim(current));
            current.clear();
        }
        else
        {
            current += c;
        }
    }

    if (!current.empty())
        result.push_back(trim(current));

    return result;
}

vector<string> extractProblemIndexCommentValues(const string& fileText)
{
    vector<string> result;

    size_t start = fileText.find("problem index corresponding to the x value");

    if (start == string::npos)
        return result;

    size_t commentEnd = fileText.find("*/", start);

    if (commentEnd == string::npos)
        return result;

    string block = fileText.substr(start, commentEnd - start);

    stringstream ss(block);
    string line;

    bool skipTitle = true;

    while (getline(ss, line))
    {
        line = trim(line);

        if (line.empty())
            continue;

        if (skipTitle)
        {
            skipTitle = false;
            continue;
        }

        result.push_back(line);
    }

    return result;
}

string extractArrayContent(const string& fileText, const string& arrayName)
{
    size_t namePos = fileText.find(arrayName);

    if (namePos == string::npos)
        return "";

    size_t start = fileText.find('{', namePos);

    if (start == string::npos)
        return "";

    int depth = 0;

    for (size_t i = start; i < fileText.size(); i++)
    {
        if (fileText[i] == '{')
            depth++;
        else if (fileText[i] == '}')
            depth--;

        if (depth == 0)
        {
            return fileText.substr(start + 1, i - start - 1);
        }
    }

    return "";
}

int extractDensitySize(const string& fileText)
{
    regex pattern(R"(double\s+density_cases\s*\[\s*\d+\s*\]\s*\[\s*(\d+)\s*\])");
    smatch match;

    if (regex_search(fileText, match, pattern))
    {
        return stoi(match[1]);
    }

    return 0;
}

vector<string> extractDensityRows(const string& densityContent)
{
    vector<string> rows;

    int depth = 0;
    size_t rowStart = string::npos;

    for (size_t i = 0; i < densityContent.size(); i++)
    {
        char c = densityContent[i];

        if (c == '{')
        {
            depth++;

            if (depth == 1)
                rowStart = i + 1;
        }
        else if (c == '}')
        {
            if (depth == 1 && rowStart != string::npos)
            {
                rows.push_back(trim(densityContent.substr(rowStart, i - rowStart)));
                rowStart = string::npos;
            }

            depth--;
        }
    }

    return rows;
}

void writeSelectedHeader(
    const string& outputFile,
    const vector<string>& xValues,
    const vector<string>& dxValues,
    const vector<string>& dsmValues,
    const vector<string>& problemIndexValues,
    const vector<string>& densityRows,
    const vector<int>& selectedIndexes,
    int densitySize
)
{
    ofstream out(outputFile);

    if (!out)
    {
        cerr << "Could not create output file: " << outputFile << endl;
        return;
    }

    int caseNumber = static_cast<int>(selectedIndexes.size());

    out << "#ifndef MAXIMUM_GENERIC_NAME_HH\n";
    out << "#define MAXIMUM_GENERIC_NAME_HH\n\n";

    out << "int Num_cases = " << caseNumber << ";\n\n";

    out << "/*\n";
    out << "problem index corresponding to the x value\n";

    for (int i = 0; i < caseNumber; i++)
    {
        int index = selectedIndexes[i];
        out << problemIndexValues[index] << " the problem is with" << indexValues[stoi(problemIndexValues[index])%4] << "\n";
    }

    out << "*/\n\n";

    out << "double x_cases[" << caseNumber << "] = {\n";
    for (int i = 0; i < caseNumber; i++)
    {
        int index = selectedIndexes[i];
        out << "    " << xValues[index];

        if (i < caseNumber - 1)
            out << ",";

        out << "\n";
    }
    out << "};\n\n";

    out << "double dsm_cases[" << caseNumber << "] = {\n";
    for (int i = 0; i < caseNumber; i++)
    {
        int index = selectedIndexes[i];
        out << "    " << dsmValues[index];

        if (i < caseNumber - 1)
            out << ",";

        out << "\n";
    }
    out << "};\n\n";

    out << "double dx_cases[" << caseNumber << "] = {\n";
    for (int i = 0; i < caseNumber; i++)
    {
        int index = selectedIndexes[i];
        out << "    " << dxValues[index];

        if (i < caseNumber - 1)
            out << ",";

        out << "\n";
    }
    out << "};\n\n";

    out << "double density_cases[" << caseNumber << "][" << densitySize << "] = {\n";
    for (int i = 0; i < caseNumber; i++)
    {
        int index = selectedIndexes[i];

        out << "    { " << densityRows[index] << " }";

        if (i < caseNumber - 1)
            out << ",";

        out << "\n";
    }
    out << "};\n\n";

    out << "#endif\n";

    cout << "Created reduced header file: " << outputFile << endl;
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "just need 1 file, and that is the UltimateGenericName.hh file\n";
        return 1;
    }
    string ultimate = argv[1];
    int mode = stoi(argv[2]);

    string inputFile =
        "./solvingProblems/" + ultimate;

    string outputFile = "MaximumGenericName.hh";

    string fileText = readWholeFile(inputFile);
    

    if (fileText.empty())
        return 1;

    string xContent = extractArrayContent(fileText, "x_cases");
    string dxContent = extractArrayContent(fileText, "dx_cases");
    string densityContent = extractArrayContent(fileText, "density_cases");
    string dsmContent = extractArrayContent(fileText, "dsm_cases");

    int densitySize = extractDensitySize(fileText);

    vector<string> xValues = splitTopLevelCommas(xContent);
    vector<string> dxValues = splitTopLevelCommas(dxContent);
    vector<string> densityRows = extractDensityRows(densityContent);
    vector<string> problemIndexValues = extractProblemIndexCommentValues(fileText);
    vector<string> dsmValues = splitTopLevelCommas(dsmContent);

    if (xValues.empty() || dxValues.empty() || dsmValues.empty() ||
    densityRows.empty() || problemIndexValues.empty())
    {
        cerr << "Failed to extract arrays from input file.\n";
        return 1;
    }
    cout << "Cases found: " << xValues.size() << endl;

    //cout << "\nWhat do you want to do?\n";
    //cout << "1. Find one specific case\n";
    //cout << "2. Find multiple specific cases\n";
    //cout << "3. Find largest dsm cases\n";
    //cout << "Choice: ";

    //int mode;
    //cin >> mode;

    vector<int> selectedIndexes;

    if (mode == 1)
    {
        int userChoice;

        cout << "Select case number, starting at 1: ";
        cin >> userChoice;

        int index = userChoice - 1;

        if (index < 0 || index >= static_cast<int>(xValues.size()))
        {
            cerr << "Invalid case number.\n";
            return 1;
        }

        selectedIndexes.push_back(index);
    }
    else if (mode == 2)
    {
        int amount;

        cout << "How many cases do you want to select? ";
        cin >> amount;

        for (int i = 0; i < amount; i++)
        {
            int userChoice;

            cout << "Enter case number #" << i + 1 << ": ";
            cin >> userChoice;

            int index = userChoice - 1;

            if (index < 0 || index >= static_cast<int>(xValues.size()))
            {
                cerr << "Invalid case number: " << userChoice << endl;
                return 1;
            }

            selectedIndexes.push_back(index);
        }
    }
    else if (mode == 3)
    {
        int largestCount = DEFAULT_LARGEST_COUNT;

        cout << "How many largest dsm cases do you want? ";
        cin >> largestCount;

        vector<pair<double, int>> dsmWithIndex;

        for (int i = 0; i < static_cast<int>(dsmValues.size()); i++)
        {
            double dsm = stod(dsmValues[i]);
            dsmWithIndex.push_back({ dsm, i });
        }

        sort(
            dsmWithIndex.begin(),
            dsmWithIndex.end(),
            [](const pair<double, int>& a, const pair<double, int>& b)
            {
                return a.first > b.first;
            }
        );

        if (largestCount > static_cast<int>(dsmWithIndex.size()))
            largestCount = static_cast<int>(dsmWithIndex.size());

        for (int i = 0; i < largestCount; i++)
        {
            selectedIndexes.push_back(dsmWithIndex[i].second);

            cout << "Selected case "
                << dsmWithIndex[i].second + 1
                << " with dsm = "
                << dsmWithIndex[i].first
                << endl;
        }
    }
    else
    {
        cerr << "Invalid option.\n";
        return 1;
    }
    

    writeSelectedHeader(
        outputFile,
        xValues,
        dxValues,
        dsmValues,
        problemIndexValues,
        densityRows,
        selectedIndexes,
        densitySize
    );
}
