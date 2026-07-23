#include "../../code/include.hh"
#include "../../run/variable.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;

void print_dummy_row(ofstream& file, int number_of_bins)
{
    for (int bin = 0; bin < number_of_bins; ++bin)
    {
        if (bin > 0)
        {
            file << ",";
        }

        file << 0;
    }

    file << endl;
}

int main(int argc, char* argv[])
{
    /*
        Usage:

        calculate_load_factors.exe collision_type N_trap eps_max output.csv

        Example:

        calculate_load_factors.exe 0 201 20 load_factors.csv
    */
    if (argc != 5)
    {
        cerr
            << "Usage:\n"
            << "    " << argv[0]
            << " collision_type N_trap eps_max output.csv\n";

        return 1;
    }

    try
    {
        const int collision_type = std::stoi(argv[1]);
        const int number_of_trapezoid_bins = std::stoi(argv[2]);
        const double epsilon_maximum = std::stod(argv[3]);
        const std::string output_filename = argv[4];

        if (collision_type != 0 && collision_type != 1)
        {
            throw std::invalid_argument(
                "collision_type must be 0 for nu-nu or 1 for nu-e."
            );
        }

        if (number_of_trapezoid_bins <= 0)
        {
            throw std::invalid_argument(
                "N_trap must be greater than zero."
            );
        }

        if (epsilon_maximum <= 0.0)
        {
            throw std::invalid_argument(
                "eps_max must be greater than zero."
            );
        }

        linspace_and_gl* epsilon_grid =
            new linspace_and_gl(
                0.0,
                epsilon_maximum,
                number_of_trapezoid_bins,
                5
            );

        const int number_of_bins = number_of_trapezoid_bins + 5;

        dep_vars* load_factors = new dep_vars(number_of_bins);

        cout
            << "Calculating load factors for "
            << number_of_bins
            << " bins..."
            << endl;

        //double results[2] = {0.0, 0.0};
        double total_estimated_load = 0.0;

    for (int bin = 0; bin < number_of_bins; ++bin)
    {
        collision_integral* integrator = nullptr;

        if (collision_type == 0)
        {
            integrator =
                new nu_nu_collision(
                    bin,
                    epsilon_grid,
                    true
                );
        }
        else
        {
            integrator =
                new nu_e_collision(
                    bin,
                    epsilon_grid,
                    true,
                    IC_TCM
                );
        }

        const double bin_load =
            integrator->estimate_load();

        load_factors->set_value(bin, bin_load);
        total_estimated_load += bin_load;

        delete integrator;

        cout
            << "\rCalculated bin "
            << bin + 1
            << " of "
            << number_of_bins
            << " | load = "
            << bin_load
            << std::flush;
    }

        cout << endl;

        if (total_estimated_load <= 0.0)
        {
            delete load_factors;
            delete epsilon_grid;

            throw std::runtime_error(
                "Every calculated load factor was zero. "
                "The header cannot balance zero-load jobs."
            );
        }

        cout
            << "Total estimated load: "
            << total_estimated_load
            << endl;

        cout << "\nLoad-factor calculation finished." << endl;

        ofstream file(output_filename);

        if (!file)
        {
            delete load_factors;
            delete epsilon_grid;

            throw std::runtime_error(
                "Could not create output CSV: " + output_filename
            );
        }

        // Row 1: dummy values.
        print_dummy_row(file, number_of_bins);

        // epsilon_grid->print_csv() produces two rows:
        // Row 2: epsilon values
        // Row 3: quadrature weights
        epsilon_grid->print_csv(file);
        file << endl;

        // Row 4: load factors
        load_factors->print_csv(file);
        file << endl;

        file.close();

        delete load_factors;
        delete epsilon_grid;
    }
    catch (const std::exception& error)
    {
        cerr
            << "calculate_load_factors stopped because of an error:\n"
            << error.what()
            << endl;

        return 1;
    }
}
