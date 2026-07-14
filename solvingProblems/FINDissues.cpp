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

*/

using namespace std;

int errorSum = 0;
int errorIndex = 0;
string indexValues [4] = {"P_0", "P_x", "P_y", "P_z"};

struct dsmRecord{
    double dsm;
    string dsmText;
    string xText;
};

void mostCommond(ofstream& out, const map<string, int>& counts, const string& label)
{
    vector<pair<string, int>> items(counts.begin(), counts.end());

    sort(items.begin(), items.end(),
        [](const pair<string, int>& a, const pair<string, int>& b)
        {
            return a.second > b.second;
        });

    out << "\nMost common " << label << ":\n";

    //int limit = min(10, static_cast<int>(items.size()));

    size_t limit = min<size_t>(10, items.size());

    for (size_t i = 0; i < limit; i++)  
    {
        if (label == "index values")
        {
            int indexNum = stoi(items[i].first);
            int modIndex = indexNum % 4;

            out << label << " = " << items[i].first
                << " appeared " << items[i].second << " times\t\t"
                << modIndex << " is the index value thus "
                << indexValues[modIndex] << " is the problem\n";
        }
        else
        {
            out << label << " = " << items[i].first
                << " appeared " << items[i].second << " times\n";
        }

        if (label == "index values")
        {
            errorSum += items[i].second;
        }
    }

}

void indexRanges(ofstream& out, const map<int, int>& indexCounts)
{
    out << "\nCommon repeated index ranges:\n";

    vector<int> indexes;

    for (const auto& pair : indexCounts)
    {
        if (pair.second > 1)
            indexes.push_back(pair.first);
    }

    if (indexes.empty())
    {
        out << "No repeated index ranges found.\n";
        return;
    }

    int start = indexes[0];
    int prev = indexes[0];

    for (size_t i = 1; i < indexes.size(); i++)
    {
        if (indexes[i] == prev + 1)
        {
            prev = indexes[i];
        }
        else
        {
            if (start == prev)
                out << start << "\n";
            else
                out << start << " - " << prev << "\n";

            start = indexes[i];
            prev = indexes[i];
        }
    }

    if (start == prev)
        out << start << "\n";
    else
        out << start << " - " << prev << "\n";
}

void writeCommonDsmRanges(ofstream& out, map<int, int>& dsmRoundedCounts)
{
    out << "\nCommon repeated DSM ranges rounded to 2 decimals:\n";

    vector<int> keys;

    for (const auto& pair : dsmRoundedCounts)
    {
        if (pair.second > 1)
            keys.push_back(pair.first);
    }

    if (keys.empty())
    {
        out << "No repeated DSM rounded ranges found.\n";
        return;
    }

    int start = keys[0];
    int prev = keys[0];

    for (size_t i = 1; i < keys.size(); i++)
    {
        if (keys[i] == prev + 1)
        {
            prev = keys[i];
        }
        else
        {
            out << fixed << setprecision(2);

            if (start == prev)
                out << start / 100.0 << "\n";
            else
                out << start / 100.0 << " - " << prev / 100.0 << "\n";

            start = keys[i];
            prev = keys[i];
        }
    }

    out << fixed << setprecision(2);

    if (start == prev)
        out << start / 100.0 << "\n";
    else
        out << start / 100.0 << " - " << prev / 100.0 << "\n";
}

void writeTopDsmValues(ofstream& out, vector<dsmRecord> records)
{
    out << "\nTop 10 largest DSM values:\n";

    if (records.empty())
    {
        out << "No DSM values found.\n";
        return;
    }

    sort(records.begin(), records.end(),
        [](const dsmRecord& a, const dsmRecord& b)
        {
            return a.dsm > b.dsm;
        });

    size_t limit = min<size_t>(10, records.size());

    for (size_t i = 0; i < limit; i++)
    {
        out << "dsm = " << records[i].dsmText
            << ", x = " << records[i].xText << "\n";
    }
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "mut have 1 input file, and that is a slurm file\n";
        return 1;
    }
    string slurm = argv[1];
    string filename = "./solvingProblems/" + slurm;
    //slurm-18329938-Copy   output data 1
    //slurm-18475761        output data 2
    //slurm-18475761.txt

    ifstream file(filename);

     if (!file)
    {
        cerr << "Could not open file: " << filename << endl;
        return 1;
    }

    ofstream out("error_info.txt");

    if (!out)
    {
        cerr << "Could not create error_info.txt\n";
        return 1;
    }

    regex indexPattern(R"(index\s*=\s*([0-9]+))");
    regex xPattern(R"(x\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?))");
    regex dxPattern(R"(dx\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?))");
    regex dsmPattern(R"(dsm\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?))");


    map<string, int> indexTextCounts;
    map<string, int> xTextCounts;
    map<string, int> dxTextCounts;

    map<string, int> dsmRoundedTextCounts;
    map<int, int> dsmRoundedNumberCounts;

    vector<dsmRecord> dsmR;


    map<int, int> indexNumberCounts;

    string line;
    string currentX = "unknown";

    while (getline(file, line))
    {
        smatch match;

        if (regex_search(line, match, xPattern))
        {
            currentX = match[1];
            xTextCounts[currentX]++;
        }

        if (regex_search(line, match, dxPattern))
        {
            string dxText = match[1];
            dxTextCounts[dxText]++;
        }

        if (regex_search(line, match, indexPattern))
        {
            string indexText = match[1];
            int indexNum = stoi(indexText);

            indexTextCounts[indexText]++;
            indexNumberCounts[indexNum]++;
        }

        if (regex_search(line, match, dsmPattern))
        {
            string dsmText = match[1];
            double dsmValue = stod(dsmText);

            double rounded = round(dsmValue * 10.0) / 10.0;
            int roundedKey = static_cast<int>(round(rounded * 10.0));

            ostringstream roundedStream;
            roundedStream << fixed << setprecision(1) << rounded;
            string roundedText = roundedStream.str();

            dsmRoundedTextCounts[roundedText]++;
            dsmRoundedNumberCounts[roundedKey]++;

            dsmR.push_back({ dsmValue, dsmText, currentX });
        }
    }
    out << "Error analysis for file: \n";
    out << "=========================================\n";

    mostCommond(out, indexTextCounts, "index values");
    writeTopDsmValues(out, dsmR);
    
    mostCommond(out, xTextCounts, "x values");
    mostCommond(out, dxTextCounts, "dx values");
    mostCommond(out, dsmRoundedTextCounts, "DSM values rounded to 1 decimals");


    indexRanges(out, indexNumberCounts);
    //mostCommond(out, dsmRoundedTextCounts, "DSM values rounded to 1 decimal");
    writeCommonDsmRanges(out, dsmRoundedNumberCounts);


    out << "\nAmount of errors found = " << errorSum << " \n";

    out << "\nAnalysis finished.\n";

    file.close();
    out.close();

    cout << "Done. Results written to error_info.txt\n";

    return 0;
}