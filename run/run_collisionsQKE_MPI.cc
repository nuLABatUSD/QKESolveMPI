#include "../code/include.hh"

#include <mpi.h>

#include "../solvingProblems/WorkerHeaderCreation/headers/collision_distribution.hh"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

namespace
{
    /*
        The number of rows is taken directly from the generated header.

        If collision_distribution.hh contains:
            int core_jobs[12][...]
        then distribution_core_count == 12.

        If it contains:
            int core_jobs[128][...]
        then distribution_core_count == 128.
    */
    const int distribution_core_count =
        sizeof(core_jobs) / sizeof(core_jobs[0]);

    int** make_core_job_rows(int numprocs)
    {
        int** rows = new int*[numprocs];

        for (int rank = 0; rank < numprocs; ++rank)
        {
            rows[rank] = core_jobs[rank];
        }

        return rows;
    }

    long long get_total_load(long long local_load)
    {
        long long total_load = 0;

        MPI_Reduce(
            &local_load,
            &total_load,
            1,
            MPI_LONG_LONG,
            MPI_SUM,
            0,
            MPI_COMM_WORLD
        );

        return total_load;
    }

    double run_compute_R(
        collisions* collision,
        std::vector<double>& results)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        const auto start =
            high_resolution_clock::now();

        collision->compute_R(
            32.0,
            32.0,
            results.data()
        );

        MPI_Barrier(MPI_COMM_WORLD);

        const auto stop =
            high_resolution_clock::now();

        return duration_cast<milliseconds>(
            stop - start
        ).count() / 1000.0;
    }
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int myid = 0;
    int numprocs = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /*
        argv[1] is NOT the number of cores.

        The shell script supplies the number of cores to mpiexec.
        argv[1] is only the test mode selected interactively by the script:

            0 = original constructor
            1 = optimized constructor
            2 = run both and compare
    */
    int mode = 0;

    if (argc >= 2)
    {
        mode = std::atoi(argv[1]);
    }

    if (mode < 0 || mode > 2)
    {
        if (myid == 0)
        {
            cerr
                << "Invalid mode.\n"
                << "0 = original\n"
                << "1 = optimized\n"
                << "2 = compare\n";
        }

        MPI_Finalize();
        return 1;
    }

    if (numprocs < 2)
    {
        if (myid == 0)
        {
            cerr
                << "At least two MPI ranks are required.\n";
        }

        MPI_Finalize();
        return 1;
    }

    /*
        The generated distribution must match the MPI core count.

        This prevents a 12-core run from reading past a 12-row header,
        and prevents a 128-core header from silently being used with
        a different number of runtime ranks.
    */
    if ((mode == 1 || mode == 2) &&
        distribution_core_count != numprocs)
    {
        if (myid == 0)
        {
            cerr
                << "Distribution/core mismatch.\n"
                << "The compiled collision_distribution.hh contains "
                << distribution_core_count
                << " core rows, but MPI started "
                << numprocs
                << " ranks.\n\n"
                << "Regenerate collision_distribution.hh for "
                << numprocs
                << " cores before running optimized mode.\n";
        }

        MPI_Finalize();
        return 1;
    }

    linspace_and_gl* eps =
        new linspace_and_gl(
            0.0,
            20.0,
            201,
            5
        );

    const int N_bins =
        eps->get_length();

    const int result_count =
        4 * N_bins;

    collisions* original_collision = nullptr;
    collisions* optimized_collision = nullptr;

    int** core_job_rows = nullptr;

    /*
        For a fair comparison, the original constructor is nu-nu only.

        The optimized constructor also creates only nu_nu_collision
        integrators.
    */
    if (mode == 0 || mode == 2)
    {
        original_collision =
            new collisions(
                myid,
                numprocs,
                eps,
                true,
                false,
                false
            );
    }

    if (mode == 1 || mode == 2)
    {
        core_job_rows =
            make_core_job_rows(numprocs);

        optimized_collision =
            new collisions(
                myid,
                numprocs,
                eps,
                core_job_rows,
                max_jobs
            );
    }

    //==============================================================
    // MODE 0: ORIGINAL CONSTRUCTOR
    //
    // Constructor/load test only.
    // No compute_R(), matching the lightweight test requested.
    //==============================================================

    if (mode == 0)
    {
        const long long local_load =
            static_cast<long long>(
                original_collision->get_load_value()
            );

        cout
            << "Rank "
            << myid
            << " original estimated load = "
            << local_load
            << endl;

        const long long total_load =
            get_total_load(local_load);

        if (myid == 0)
        {
            cout
                << "\nMode: original collision distribution\n"
                << "MPI cores: "
                << numprocs
                << "\n"
                << "Epsilon bins: "
                << N_bins
                << "\n"
                << "Total estimated workload: "
                << total_load
                << "\n";
        }
    }

    //==============================================================
    // MODE 1: OPTIMIZED CONSTRUCTOR
    //
    // Constructor/load test only.
    // Uses core_jobs directly from collision_distribution.hh.
    // No compute_R().
    //==============================================================

    else if (mode == 1)
    {
        const long long local_load =
            static_cast<long long>(
                optimized_collision->get_load_value()
            );

        cout
            << "Rank "
            << myid
            << " optimized estimated load = "
            << local_load
            << endl;

        const long long total_load =
            get_total_load(local_load);

        if (myid == 0)
        {
            cout
                << "\nMode: optimized collision distribution\n"
                << "MPI cores: "
                << numprocs
                << "\n"
                << "Distribution header cores: "
                << distribution_core_count
                << "\n"
                << "Epsilon bins: "
                << N_bins
                << "\n"
                << "Maximum jobs per worker: "
                << max_jobs
                << "\n"
                << "Total estimated workload: "
                << total_load
                << "\n";
        }
    }

    //==============================================================
    // MODE 2: VERIFY BOTH
    //
    // This is the optional heavy test.
    // It compares total workload and all compute_R() results.
    //==============================================================

    else
    {
        const long long original_local_load =
            static_cast<long long>(
                original_collision->get_load_value()
            );

        const long long optimized_local_load =
            static_cast<long long>(
                optimized_collision->get_load_value()
            );

        cout
            << "Rank "
            << myid
            << ": original load = "
            << original_local_load
            << ", optimized load = "
            << optimized_local_load
            << endl;

        const long long original_total_load =
            get_total_load(original_local_load);

        const long long optimized_total_load =
            get_total_load(optimized_local_load);

        std::vector<double> original_results(
            result_count,
            0.0
        );

        std::vector<double> optimized_results(
            result_count,
            0.0
        );

        const double original_time =
            run_compute_R(
                original_collision,
                original_results
            );

        const double optimized_time =
            run_compute_R(
                optimized_collision,
                optimized_results
            );

        if (myid == 0)
        {
            const double absolute_tolerance = 1.0e-18;
            const double relative_tolerance = 1.0e-10;

            int exact_matches = 0;
            int within_tolerance = 0;
            int failed_values = 0;

            double maximum_absolute_difference = 0.0;
            double maximum_relative_difference = 0.0;

            int max_difference_index = -1;

            for (int index = 0;
                 index < result_count;
                 ++index)
            {
                const double original_value =
                    original_results[index];

                const double optimized_value =
                    optimized_results[index];

                if (original_value == optimized_value)
                {
                    ++exact_matches;
                }

                const double absolute_difference =
                    std::abs(
                        original_value -
                        optimized_value
                    );

                const double scale =
                    std::max(
                        std::abs(original_value),
                        std::abs(optimized_value)
                    );

                const double relative_difference =
                    scale > 0.0
                        ? absolute_difference / scale
                        : 0.0;

                const bool matches =
                    absolute_difference <= absolute_tolerance ||
                    absolute_difference <=
                        relative_tolerance * scale;

                if (matches)
                {
                    ++within_tolerance;
                }
                else
                {
                    ++failed_values;
                }

                if (absolute_difference >
                    maximum_absolute_difference)
                {
                    maximum_absolute_difference =
                        absolute_difference;

                    max_difference_index =
                        index;
                }

                if (relative_difference >
                    maximum_relative_difference)
                {
                    maximum_relative_difference =
                        relative_difference;
                }
            }

            const bool workload_matches =
                original_total_load ==
                optimized_total_load;

            const bool results_match =
                failed_values == 0;

            cout << std::setprecision(17);

            cout
                << "\n================ Comparison ================\n"
                << "MPI cores:                    "
                << numprocs
                << "\n"
                << "Distribution header cores:    "
                << distribution_core_count
                << "\n"
                << "Original total workload:      "
                << original_total_load
                << "\n"
                << "Optimized total workload:     "
                << optimized_total_load
                << "\n"
                << "Workload difference:          "
                << original_total_load -
                   optimized_total_load
                << "\n"
                << "Original compute_R time:      "
                << original_time
                << " seconds\n"
                << "Optimized compute_R time:     "
                << optimized_time
                << " seconds\n";

            if (optimized_time > 0.0)
            {
                cout
                    << "Measured speedup:             "
                    << original_time /
                       optimized_time
                    << "x\n";
            }

            cout
                << "\nR values compared:            "
                << result_count
                << "\n"
                << "Exact matches:                "
                << exact_matches
                << "\n"
                << "Within tolerance:             "
                << within_tolerance
                << "\n"
                << "Outside tolerance:            "
                << failed_values
                << "\n"
                << "Maximum absolute difference:  "
                << maximum_absolute_difference
                << "\n"
                << "Maximum relative difference:  "
                << maximum_relative_difference
                << "\n";

            if (max_difference_index >= 0)
            {
                const int block =
                    max_difference_index / N_bins;

                const int bin =
                    max_difference_index % N_bins;

                cout
                    << "Largest difference R block:   "
                    << block
                    << "\n"
                    << "Largest difference bin:       "
                    << bin
                    << "\n"
                    << "Largest difference epsilon:   "
                    << eps->get_value(bin)
                    << "\n";
            }

            cout
                << "\nWorkload check: "
                << (workload_matches ? "PASS" : "FAIL")
                << "\n"
                << "R-value check:  "
                << (results_match ? "PASS" : "FAIL")
                << "\n"
                << "Overall check:  "
                << (
                    workload_matches &&
                    results_match
                        ? "PASS"
                        : "FAIL"
                )
                << "\n"
                << "============================================\n";
        }
    }

    delete optimized_collision;
    delete original_collision;

    delete[] core_job_rows;

    delete eps;

    MPI_Finalize();

    return 0;
}
