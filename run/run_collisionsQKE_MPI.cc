#include "../code/include.hh"

#include <mpi.h>

#include "../solvingProblems/WorkerHeaderCreation/headers/collision_distribution_128_cores_201_trap.hh"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

int main1(int argc, char* argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        
    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);
    int N_bins = eps->get_length();
    

    collisions* C_MPI = new collisions(myid, numprocs, eps, true, true, false);
    //int[] [] worker = collision_distribution_128_cores_201_trap.core_jobs;
    //int** worker = core_jobs;
    //collisions* C_MPI = new collisions(myid, numprocs, eps, core_jobs, 128);
        
    double* R_values = new double[N_bins * 4];

    auto start = high_resolution_clock::now();
    
    C_MPI->compute_R(32., 32., R_values); 
        
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    cout << myid << ", " << C_MPI->get_load_value() << endl;
    
    if(myid == 0){
        for(int i = 0; i < eps->get_length(); i++)
            cout << eps->get_value(i) << ", " << R_values[i] << ", " << R_values[2*206+i] << endl;
        cout << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;
    }

    delete[] R_values;
    delete C_MPI;
        
    delete eps;
    
    MPI_Finalize();
    return 0;
}

namespace
{
    constexpr int source_core_count = 128;

    struct RuntimeDistribution
    {
        int** jobs = nullptr;
        int max_jobs_per_worker = 0;
        int recovered_jobs = 0;
    };

    struct WorkloadSummary
    {
        long long total = 0;
        long long minimum = 0;
        long long maximum = 0;
    };

    struct ResultSummary
    {
        std::size_t exact_matches = 0;
        std::size_t within_tolerance = 0;
        std::size_t outside_tolerance = 0;
        std::size_t non_finite = 0;

        double maximum_absolute_difference = 0.0;
        double maximum_relative_difference = 0.0;

        int maximum_absolute_index = -1;
        int maximum_relative_index = -1;
    };

    bool nearly_equal(
        double first,
        double second,
        double absolute_tolerance,
        double relative_tolerance)
    {
        if (!std::isfinite(first) ||
            !std::isfinite(second))
        {
            return false;
        }

        const double absolute_difference =
            std::abs(first - second);

        if (absolute_difference <= absolute_tolerance)
        {
            return true;
        }

        const double scale =
            std::max(
                std::abs(first),
                std::abs(second));

        return absolute_difference <=
               relative_tolerance * scale;
    }

    RuntimeDistribution create_runtime_distribution(
        int numprocs,
        int number_of_bins,
        int myid)
    {
        RuntimeDistribution distribution;

        const int runtime_worker_count =
            numprocs - 1;

        std::vector<std::vector<int>>
            worker_jobs(numprocs);

        std::vector<int> all_jobs;

        for (int source_rank = 1;
             source_rank < source_core_count;
             ++source_rank)
        {
            const int runtime_rank =
                1 + (
                    (source_rank - 1) %
                    runtime_worker_count);

            for (int slot = 0;
                 slot < max_jobs;
                 ++slot)
            {
                const int job =
                    core_jobs[source_rank][slot];

                if (job < 0)
                {
                    continue;
                }

                worker_jobs[runtime_rank]
                    .push_back(job);

                all_jobs.push_back(job);
            }
        }

        std::vector<int> sorted_jobs =
            all_jobs;

        std::sort(
            sorted_jobs.begin(),
            sorted_jobs.end());

        const int expected_jobs =
            2 * number_of_bins;

        const auto duplicate =
            std::adjacent_find(
                sorted_jobs.begin(),
                sorted_jobs.end());

        const int minimum_job =
            sorted_jobs.empty()
                ? -1
                : sorted_jobs.front();

        const int maximum_job =
            sorted_jobs.empty()
                ? -1
                : sorted_jobs.back();

        if (duplicate != sorted_jobs.end() ||
            static_cast<int>(
                sorted_jobs.size()) != expected_jobs ||
            minimum_job != 0 ||
            maximum_job != expected_jobs - 1)
        {
            if (myid == 0)
            {
                cerr
                    << "Invalid optimized job distribution.\n"
                    << "Expected "
                    << expected_jobs
                    << " unique jobs numbered 0 through "
                    << expected_jobs - 1
                    << ".\n"
                    << "Recovered: "
                    << sorted_jobs.size()
                    << "\nMinimum: "
                    << minimum_job
                    << "\nMaximum: "
                    << maximum_job
                    << '\n';

                if (duplicate != sorted_jobs.end())
                {
                    cerr
                        << "Duplicate job: "
                        << *duplicate
                        << '\n';
                }
            }

            return distribution;
        }

        for (int rank = 1;
             rank < numprocs;
             ++rank)
        {
            distribution.max_jobs_per_worker =
                std::max(
                    distribution.max_jobs_per_worker,
                    static_cast<int>(
                        worker_jobs[rank].size()));
        }

        distribution.jobs =
            new int*[numprocs];

        for (int rank = 0;
             rank < numprocs;
             ++rank)
        {
            distribution.jobs[rank] =
                new int[
                    distribution.max_jobs_per_worker];

            std::fill(
                distribution.jobs[rank],
                distribution.jobs[rank] +
                    distribution.max_jobs_per_worker,
                -1);
        }

        for (int rank = 1;
             rank < numprocs;
             ++rank)
        {
            for (std::size_t slot = 0;
                 slot < worker_jobs[rank].size();
                 ++slot)
            {
                distribution.jobs[rank][slot] =
                    worker_jobs[rank][slot];
            }
        }

        distribution.recovered_jobs =
            static_cast<int>(all_jobs.size());

        return distribution;
    }

    void destroy_runtime_distribution(
        RuntimeDistribution& distribution,
        int numprocs)
    {
        if (distribution.jobs == nullptr)
        {
            return;
        }

        for (int rank = 0;
             rank < numprocs;
             ++rank)
        {
            delete[] distribution.jobs[rank];
        }

        delete[] distribution.jobs;
        distribution.jobs = nullptr;
    }

    WorkloadSummary collect_workload_summary(
        long long local_load,
        int myid)
    {
        WorkloadSummary summary;

        const long long minimum_candidate =
            myid == 0
                ? LLONG_MAX
                : local_load;

        MPI_Reduce(
            &local_load, &summary.total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(
            &minimum_candidate, &summary.minimum, 1, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);

        MPI_Reduce(
            &local_load, &summary.maximum, 1, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

        return summary;
    }

    ResultSummary compare_results(
        const std::vector<double>& original,
        const std::vector<double>& optimized,
        double absolute_tolerance,
        double relative_tolerance)
    {
        ResultSummary summary;

        for (std::size_t index = 0;
             index < original.size();
             ++index)
        {
            const double original_value =
                original[index];

            const double optimized_value =
                optimized[index];

            if (!std::isfinite(original_value) ||
                !std::isfinite(optimized_value))
            {
                ++summary.non_finite;
                ++summary.outside_tolerance;
                continue;
            }

            if (original_value == optimized_value)
            {
                ++summary.exact_matches;
            }

            const double absolute_difference =
                std::abs(
                    original_value -
                    optimized_value);

            const double scale =
                std::max(
                    std::abs(original_value),
                    std::abs(optimized_value));

            const double relative_difference =
                scale > 0.0
                    ? absolute_difference / scale
                    : 0.0;

            if (absolute_difference >
                summary.maximum_absolute_difference)
            {
                summary.maximum_absolute_difference =
                    absolute_difference;

                summary.maximum_absolute_index =
                    static_cast<int>(index);
            }

            if (relative_difference >
                summary.maximum_relative_difference)
            {
                summary.maximum_relative_difference =
                    relative_difference;

                summary.maximum_relative_index =
                    static_cast<int>(index);
            }

            if (nearly_equal(
                    original_value,
                    optimized_value,
                    absolute_tolerance,
                    relative_tolerance))
            {
                ++summary.within_tolerance;
            }
            else
            {
                ++summary.outside_tolerance;
            }
        }

        return summary;
    }

    double run_collision(
        collisions* collision,
        std::vector<double>& results)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        const auto start =
            high_resolution_clock::now();

        collision->compute_R(
            32.0,
            32.0,
            results.data());

        MPI_Barrier(MPI_COMM_WORLD);

        const auto stop =
            high_resolution_clock::now();

        return duration_cast<milliseconds>(
                   stop - start).count() /
               1000.0;
    }

    void print_selected_results(
        linspace_and_gl* eps,
        const std::vector<double>& results,
        int number_of_bins)
    {
        for (int index = 0;
             index < number_of_bins;
             ++index)
        {
            cout
                << eps->get_value(index)
                << ", "
                << results[index]
                << ", "
                << results[
                    (2 * number_of_bins) + index]
                << '\n';
        }
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
    ------------------------------------------------------
     Select operating mode
    
     0 = original nu-nu constructor
     1 = optimized constructor
     2 = run both and compare workload and R values
    ------------------------------------------------------
    */

    int process_mode = 2;

    if (argc >= 2)
    {
        process_mode =
            std::stoi(argv[1]);
    }

    if (process_mode < 0 ||
        process_mode > 2)
    {
        if (myid == 0)
        {
            cerr
                << "Usage:\n"
                << "    program.exe 0\n"
                << "        Run original nu-nu distribution.\n"
                << "    program.exe 1\n"
                << "        Run optimized distribution.\n"
                << "    program.exe 2 [absolute_tolerance] "
                << "[relative_tolerance]\n"
                << "        Run both and compare results.\n";
        }

        MPI_Finalize();
        return 1;
    }

    double absolute_tolerance =
        1.0e-18;

    double relative_tolerance =
        1.0e-10;

    try
    {
        if (argc >= 3)
        {
            absolute_tolerance =
                std::stod(argv[2]);
        }

        if (argc >= 4)
        {
            relative_tolerance =
                std::stod(argv[3]);
        }
    }
    catch (const std::exception& error)
    {
        if (myid == 0)
        {
            cerr
                << "Invalid tolerance: "
                << error.what()
                << '\n';
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

    linspace_and_gl* eps =
        new linspace_and_gl(0.0, 20.0, 201, 5);

    const int number_of_bins =
        eps->get_length();

    const int result_count =
        4 * number_of_bins;

    RuntimeDistribution runtime_distribution;

    if (process_mode == 1 ||
        process_mode == 2)
    {
        runtime_distribution =
            create_runtime_distribution(
                numprocs,
                number_of_bins,
                myid);

        if (runtime_distribution.jobs == nullptr)
        {
            delete eps;
            MPI_Finalize();
            return 1;
        }
    }

    collisions* original_collision =
        nullptr;

    collisions* optimized_collision =
        nullptr;

    if (process_mode == 0 ||
        process_mode == 2)
    {
        original_collision =
            new collisions(
                myid, numprocs, eps, true, false, false);
    }

    if (process_mode == 1 ||
        process_mode == 2)
    {
        optimized_collision =
            new collisions(
                myid, numprocs, eps, runtime_distribution.jobs, runtime_distribution.max_jobs_per_worker);
    }

    if (process_mode == 0)
    {
        std::vector<double> results(
            result_count,
            0.0);

        const long long local_load =
            static_cast<long long>(
                original_collision
                    ->get_load_value());

        cout
            << "Rank "
            << myid
            << " original estimated load = "
            << local_load
            << '\n';

        const double elapsed_seconds =
            run_collision(
                original_collision,
                results);

        if (myid == 0)
        {
            cout
                << "\nDistribution: original nu-nu\n"
                << "Total cores: "
                << numprocs
                << '\n'
                << "Number of epsilon bins: "
                << number_of_bins
                << '\n'
                << "Time elapsed: "
                << elapsed_seconds
                << " seconds\n\n";

            print_selected_results(
                eps,
                results,
                number_of_bins);
        }
    }
    else if (process_mode == 1)
    {
        std::vector<double> results(
            result_count,
            0.0);

        const long long local_load =
            static_cast<long long>(
                optimized_collision
                    ->get_load_value());

        cout
            << "Rank "
            << myid
            << " optimized estimated load = "
            << local_load
            << '\n';

        const double elapsed_seconds =
            run_collision(
                optimized_collision,
                results);

        if (myid == 0)
        {
            cout
                << "\nDistribution: optimized\n"
                << "Static source ranks: "
                << source_core_count
                << '\n'
                << "Runtime cores: "
                << numprocs
                << '\n'
                << "Recovered jobs: "
                << runtime_distribution
                    .recovered_jobs
                << '\n'
                << "Runtime maximum jobs per worker: "
                << runtime_distribution
                    .max_jobs_per_worker
                << '\n'
                << "Time elapsed: "
                << elapsed_seconds
                << " seconds\n\n";

            print_selected_results(
                eps,
                results,
                number_of_bins);
        }
    }
    else
    {
        std::vector<double> original_results(
            result_count,
            0.0);

        std::vector<double> optimized_results(
            result_count,
            0.0);

        const long long original_local_load =
            static_cast<long long>(
                original_collision
                    ->get_load_value());

        const long long optimized_local_load =
            static_cast<long long>(
                optimized_collision
                    ->get_load_value());

        cout
            << "Rank "
            << myid
            << ": original load = "
            << original_local_load
            << ", optimized load = "
            << optimized_local_load
            << '\n';

        const WorkloadSummary original_workload =
            collect_workload_summary(
                original_local_load,
                myid);

        const WorkloadSummary optimized_workload =
            collect_workload_summary(
                optimized_local_load,
                myid);

        const double original_seconds =
            run_collision(
                original_collision,
                original_results);

        const double optimized_seconds =
            run_collision(
                optimized_collision,
                optimized_results);

        if (myid == 0)
        {
            const ResultSummary comparison =
                compare_results(
                    original_results,
                    optimized_results,
                    absolute_tolerance,
                    relative_tolerance);

            const long long workload_difference =
                original_workload.total -
                optimized_workload.total;

            cout << std::setprecision(17);

            cout
                << "\n================ Verification ================\n"
                << "Runtime MPI ranks:           "
                << numprocs
                << '\n'
                << "Recovered optimized jobs:    "
                << runtime_distribution
                    .recovered_jobs
                << '\n'
                << "Original total workload:     "
                << original_workload.total
                << '\n'
                << "Optimized total workload:    "
                << optimized_workload.total
                << '\n'
                << "Workload difference:         "
                << workload_difference
                << '\n'
                << "Original maximum workload:   "
                << original_workload.maximum
                << '\n'
                << "Optimized maximum workload:  "
                << optimized_workload.maximum
                << '\n'
                << "Original elapsed time:       "
                << original_seconds
                << " seconds\n"
                << "Optimized elapsed time:      "
                << optimized_seconds
                << " seconds\n";

            if (optimized_seconds > 0.0)
            {
                cout
                    << "Measured speedup:            "
                    << original_seconds /
                       optimized_seconds
                    << "x\n";
            }

            cout
                << "\nR values compared:           "
                << original_results.size()
                << '\n'
                << "Exact matches:               "
                << comparison.exact_matches
                << '\n'
                << "Within tolerance:            "
                << comparison.within_tolerance
                << '\n'
                << "Outside tolerance:           "
                << comparison.outside_tolerance
                << '\n'
                << "Non-finite comparisons:      "
                << comparison.non_finite
                << '\n'
                << "Maximum absolute difference: "
                << comparison
                    .maximum_absolute_difference
                << '\n'
                << "Maximum relative difference: "
                << comparison
                    .maximum_relative_difference
                << '\n';

            const bool workload_pass =
                original_workload.total ==
                optimized_workload.total;

            const bool result_pass =
                comparison.outside_tolerance == 0;

            cout
                << "\nWorkload result: "
                << (
                    workload_pass
                        ? "PASS"
                        : "FAIL")
                << '\n'
                << "R-value result:  "
                << (
                    result_pass
                        ? "PASS"
                        : "FAIL")
                << '\n'
                << "Overall result:  "
                << (
                    workload_pass &&
                    result_pass
                        ? "PASS"
                        : "FAIL")
                << "\n==============================================\n";
        }
    }

    delete optimized_collision;
    delete original_collision;

    destroy_runtime_distribution(
        runtime_distribution,
        numprocs);

    delete eps;

    MPI_Finalize();

    return 0;
}
