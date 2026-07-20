#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <filesystem>

#include <algorithm>
#include <functional>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>

#include "../../code/include.hh"
#include "../../run/calculate_diagnostics.cc"


// DATA TYPES

// Use long long so combined worker loads are not limited to the
// smaller range of a 32-bit int.
using LoadValue = long long;


struct Task
{
    // Collision-integral identifier.
    int job_id = -1;

    // Physical epsilon-bin index.
    int epsilon_bin = -1;

    // Estimated computational load.
    LoadValue load = 0;
};

// assignments[worker_index] contains all jobs for that worker.
using Assignments = std::vector<std::vector<Task>>;

enum class DistributionMethod
{
    Wall,
    Snake,
    Mountain,
    LargestFirst,
    HighLowSnake
};

struct DistributionResult
{
    DistributionMethod method = DistributionMethod::Wall;

    Assignments assignments;
    std::vector<LoadValue> worker_loads;

    LoadValue maximum_load = 0;
    LoadValue minimum_load = 0;

    double average_load = 0.0;
    double efficiency = 0.0;
};


// VALIDATION AND GENERAL HELPERS

void validateDistributionInputs(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    if (worker_count <= 0)
    {
        throw std::invalid_argument(
            "worker_count must be greater than zero."
        );
    }

    if (tasks.empty())
    {
        throw std::invalid_argument(
            "The task list cannot be empty."
        );
    }

    for (const Task& task : tasks)
    {
        if (task.load < 0)
        {
            throw std::invalid_argument(
                "Task loads cannot be negative."
            );
        }
    }
}

Assignments createEmptyAssignments(int worker_count)
{
    if (worker_count <= 0)
    {
        throw std::invalid_argument(
            "worker_count must be greater than zero."
        );
    }

    return Assignments(
        static_cast<std::size_t>(worker_count)
    );
}

std::vector<LoadValue> calculateWorkerLoads(
    const Assignments& assignments
)
{
    std::vector<LoadValue> worker_loads;
    worker_loads.reserve(assignments.size());

    for (const std::vector<Task>& worker_tasks : assignments)
    {
        LoadValue total_load = 0;

        for (const Task& task : worker_tasks)
        {
            total_load += task.load;
        }

        worker_loads.push_back(total_load);
    }

    return worker_loads;
}


DistributionResult calculateDistributionResult(
    DistributionMethod method,
    const Assignments& assignments
)
{
    DistributionResult result;

    result.method = method;
    result.assignments = assignments;
    result.worker_loads = calculateWorkerLoads(assignments);

    if (result.worker_loads.empty())
    {
        return result;
    }

    const auto minimum_and_maximum = std::minmax_element(
        result.worker_loads.begin(),
        result.worker_loads.end()
    );

    result.minimum_load = *minimum_and_maximum.first;
    result.maximum_load = *minimum_and_maximum.second;

    const LoadValue total_load = std::accumulate(
        result.worker_loads.begin(),
        result.worker_loads.end(),
        LoadValue{0}
    );

    result.average_load =
        static_cast<double>(total_load)
        / static_cast<double>(result.worker_loads.size());

    if (result.maximum_load > 0)
    {
        result.efficiency =
            result.average_load
            / static_cast<double>(result.maximum_load);
    }

    return result;
}

std::string distributionMethodName(
    DistributionMethod method
)
{
    switch (method)
    {
        case DistributionMethod::Wall:
            return "Wall";

        case DistributionMethod::Snake:
            return "Snake";

        case DistributionMethod::Mountain:
            return "Mountain";

        case DistributionMethod::LargestFirst:
            return "Largest-first";

        case DistributionMethod::HighLowSnake:
            return "High-low snake";
    }

    return "Unknown";
}


// METHOD 1: WALL


Assignments distributeWall(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    validateDistributionInputs(tasks, worker_count);

    Assignments assignments =
        createEmptyAssignments(worker_count);

    for (
        std::size_t task_index = 0;
        task_index < tasks.size();
        ++task_index
    )
    {
        const int worker_index =
            static_cast<int>(task_index)
            % worker_count;

        assignments[worker_index].push_back(
            tasks[task_index]
        );
    }

    return assignments;
}


// METHOD 2: SNAKE

Assignments distributeSnake(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    validateDistributionInputs(tasks, worker_count);

    Assignments assignments =
        createEmptyAssignments(worker_count);

    for (
        std::size_t task_index = 0;
        task_index < tasks.size();
        ++task_index
    )
    {
        const int integer_task_index =
            static_cast<int>(task_index);

        const int row_number =
            integer_task_index / worker_count;

        const int position_in_row =
            integer_task_index % worker_count;

        int worker_index = 0;

        if (row_number % 2 == 0)
        {
            worker_index = position_in_row;
        }
        else
        {
            worker_index =
                worker_count - 1 - position_in_row;
        }

        assignments[worker_index].push_back(
            tasks[task_index]
        );
    }

    return assignments;
}


// METHOD 3: MOUNTAIN
// CONTIGUOUS MINIMAX PARTITION

int numberOfGroupsNeeded(
    const std::vector<Task>& tasks,
    LoadValue maximum_group_load
)
{
    int group_count = 1;
    LoadValue current_load = 0;

    for (const Task& task : tasks)
    {
        if (task.load > maximum_group_load)
        {
            return static_cast<int>(tasks.size()) + 1;
        }

        // Written this way to avoid overflow from:
        //
        //     current_load + task.load
        //
        if (
            current_load
            <= maximum_group_load - task.load
        )
        {
            current_load += task.load;
        }
        else
        {
            ++group_count;
            current_load = task.load;
        }
    }

    return group_count;
}


LoadValue findContiguousLoadLimit(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    LoadValue lower_bound = 0;
    LoadValue upper_bound = 0;

    for (const Task& task : tasks)
    {
        lower_bound = std::max(
            lower_bound,
            task.load
        );

        upper_bound += task.load;
    }

    // Integer binary search for the smallest maximum group load
    // that allows the ordered task list to fit into no more than
    // worker_count contiguous groups.
    while (lower_bound < upper_bound)
    {
        const LoadValue middle =
            lower_bound
            + (upper_bound - lower_bound) / 2;

        const int groups_needed =
            numberOfGroupsNeeded(
                tasks,
                middle
            );

        if (groups_needed <= worker_count)
        {
            upper_bound = middle;
        }
        else
        {
            lower_bound = middle + 1;
        }
    }

    return lower_bound;
}


Assignments buildExactContiguousGroups(
    const std::vector<Task>& tasks,
    int worker_count,
    LoadValue load_limit
)
{
    if (
        worker_count
        > static_cast<int>(tasks.size())
    )
    {
        throw std::invalid_argument(
            "Mountain cannot create more nonempty groups "
            "than there are tasks."
        );
    }

    Assignments assignments;
    assignments.reserve(
        static_cast<std::size_t>(worker_count)
    );

    std::vector<Task> current_group;
    LoadValue current_load = 0;

    for (
        std::size_t task_index = 0;
        task_index < tasks.size();
        ++task_index
    )
    {
        const Task& task = tasks[task_index];

        const int tasks_remaining_including_current =
            static_cast<int>(
                tasks.size() - task_index
            );

        const int workers_remaining_after_current_group =
            worker_count
            - static_cast<int>(assignments.size())
            - 1;

        // If the remaining number of tasks equals the number of
        // remaining worker groups, close the current group before
        // this task so every remaining worker receives one task.
        const bool must_start_new_group =
            !current_group.empty()
            && tasks_remaining_including_current
                == workers_remaining_after_current_group;

        const bool exceeds_limit =
            !current_group.empty()
            && current_load
                > load_limit - task.load;

        const bool can_start_new_group =
            static_cast<int>(assignments.size())
            < worker_count - 1;

        if (
            can_start_new_group
            && (
                must_start_new_group
                || exceeds_limit
            )
        )
        {
            assignments.push_back(
                current_group
            );

            current_group.clear();
            current_load = 0;
        }

        current_group.push_back(task);
        current_load += task.load;
    }

    if (!current_group.empty())
    {
        assignments.push_back(current_group);
    }

    if (
        static_cast<int>(assignments.size())
        != worker_count
    )
    {
        throw std::runtime_error(
            "Mountain did not produce exactly the requested "
            "number of worker groups."
        );
    }

    return assignments;
}


Assignments distributeMountain(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    validateDistributionInputs(tasks, worker_count);

    if (
        worker_count
        > static_cast<int>(tasks.size())
    )
    {
        throw std::invalid_argument(
            "Mountain requires worker_count <= task count."
        );
    }

    const LoadValue load_limit =
        findContiguousLoadLimit(
            tasks,
            worker_count
        );

    return buildExactContiguousGroups(
        tasks,
        worker_count,
        load_limit
    );
}


// METHOD 4: LARGEST-FIRST
// LARGEST-PROCESSING-TIME GREEDY

Assignments distributeLargestFirst(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    validateDistributionInputs(tasks, worker_count);

    Assignments assignments =
        createEmptyAssignments(worker_count);

    // Each entry stores:
    //
    //     current worker load
    //     worker index
    //
    // std::greater makes this a min-heap.
    using HeapEntry = std::pair<LoadValue, int>;

    std::priority_queue<
        HeapEntry,
        std::vector<HeapEntry>,
        std::greater<HeapEntry>
    > worker_heap;

    for (
        int worker_index = 0;
        worker_index < worker_count;
        ++worker_index
    )
    {
        worker_heap.push({
            0,
            worker_index
        });
    }

    std::vector<Task> sorted_tasks = tasks;

    std::stable_sort(
        sorted_tasks.begin(),
        sorted_tasks.end(),
        [](
            const Task& left,
            const Task& right
        )
        {
            return left.load > right.load;
        }
    );

    for (const Task& task : sorted_tasks)
    {
        const HeapEntry smallest_worker =
            worker_heap.top();

        worker_heap.pop();

        const LoadValue current_load =
            smallest_worker.first;

        const int worker_index =
            smallest_worker.second;

        assignments[worker_index].push_back(task);

        worker_heap.push({
            current_load + task.load,
            worker_index
        });
    }

    return assignments;
}


// METHOD 5: HIGH-LOW SNAKE

Assignments distributeHighLowSnake(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    validateDistributionInputs(tasks, worker_count);

    std::vector<Task> sorted_tasks = tasks;

    std::stable_sort(
        sorted_tasks.begin(),
        sorted_tasks.end(),
        [](
            const Task& left,
            const Task& right
        )
        {
            return left.load < right.load;
        }
    );

    std::vector<Task> reordered_tasks;
    reordered_tasks.reserve(sorted_tasks.size());

    std::size_t left = 0;
    std::size_t right = sorted_tasks.size() - 1;

    while (left <= right)
    {
        // Largest remaining task.
        reordered_tasks.push_back(
            sorted_tasks[right]
        );

        if (right == 0)
        {
            break;
        }

        --right;

        // Smallest remaining task.
        if (left <= right)
        {
            reordered_tasks.push_back(
                sorted_tasks[left]
            );

            ++left;
        }
    }

    return distributeSnake(
        reordered_tasks,
        worker_count
    );
}


// METHOD DISPATCH

Assignments distributeTasks(
    const std::vector<Task>& tasks,
    int worker_count,
    DistributionMethod method
)
{
    switch (method)
    {
        case DistributionMethod::Wall:
            return distributeWall(
                tasks,
                worker_count
            );

        case DistributionMethod::Snake:
            return distributeSnake(
                tasks,
                worker_count
            );

        case DistributionMethod::Mountain:
            return distributeMountain(
                tasks,
                worker_count
            );

        case DistributionMethod::LargestFirst:
            return distributeLargestFirst(
                tasks,
                worker_count
            );

        case DistributionMethod::HighLowSnake:
            return distributeHighLowSnake(
                tasks,
                worker_count
            );
    }

    throw std::invalid_argument(
        "Unknown distribution method."
    );
}
/*
before you ask, switch statements are better than if. 
    if still cheks line by line even if the statement is false
    switch will just skip that entire portion of the code all together, faster
    but has more restrictions in terms of its inputs
*/


// RUN AND COMPARE DISTRIBUTION METHODS

std::vector<DistributionResult> runAllDistributionMethods(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    const std::vector<DistributionMethod> methods = {
        DistributionMethod::Wall,
        DistributionMethod::Snake,
        DistributionMethod::Mountain,
        DistributionMethod::LargestFirst,
        DistributionMethod::HighLowSnake
    };

    std::vector<DistributionResult> results;
    results.reserve(methods.size());

    for (DistributionMethod method : methods)
    {
        const Assignments assignments =
            distributeTasks(
                tasks,
                worker_count,
                method
            );

        results.push_back(
            calculateDistributionResult(
                method,
                assignments
            )
        );
    }

    return results;
}


DistributionResult chooseBestDistribution(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    std::vector<DistributionResult> results =
        runAllDistributionMethods(
            tasks,
            worker_count
        );

    const auto best_result = std::min_element(
        results.begin(),
        results.end(),
        [](
            const DistributionResult& left,
            const DistributionResult& right
        )
        {
            return left.maximum_load
                < right.maximum_load;
        }
    );

    if (best_result == results.end())
    {
        throw std::runtime_error(
            "No distribution result was created."
        );
    }

    return *best_result;
}


DistributionResult chooseBestPrimaryDistribution(
    const std::vector<Task>& tasks,
    int worker_count
)
{
    const DistributionResult mountain_result =
        calculateDistributionResult(
            DistributionMethod::Mountain,
            distributeMountain(
                tasks,
                worker_count
            )
        );

    const DistributionResult largest_first_result =
        calculateDistributionResult(
            DistributionMethod::LargestFirst,
            distributeLargestFirst(
                tasks,
                worker_count
            )
        );

    if (
        largest_first_result.maximum_load
        < mountain_result.maximum_load
    )
    {
        return largest_first_result;
    }

    // Mountain wins ties.
    return mountain_result;
}


// METHOD SELECTION

DistributionMethod chooseDistributionMethod(int worker_count)
{
    // Benchmark result supplied for this project:
    // workers 122-127 use Mountain; all others use Largest-first.
    if (worker_count >= 122 && worker_count <= 127)
    {
        return DistributionMethod::Mountain;
    }

    return DistributionMethod::LargestFirst;
}


// HEADER CREATION HELPERS

int calculateMaximumJobs(const Assignments& assignments)
{
    int maximum_jobs = 0;

    for (const std::vector<Task>& worker_tasks : assignments)
    {
        maximum_jobs = std::max(
            maximum_jobs,
            static_cast<int>(worker_tasks.size())
        );
    }

    return maximum_jobs;
}


void writeCoreJobsMatrix1(
    std::ofstream& output,
    int total_core_count,
    int maximum_jobs,
    const Assignments& assignments
)
{
    output
        << "inline constexpr int core_jobs["
        << total_core_count
        << "]["
        << maximum_jobs
        << "] = {\n";

    // MPI core 0 is the coordinator and receives no jobs.
    output << "    {";

    for (int column = 0; column < maximum_jobs; ++column)
    {
        if (column > 0)
        {
            output << ", ";
        }

        output << -1;
    }

    output << "},\n";

    // assignments[0] belongs to MPI core 1.
    for (int core_index = 1; core_index < total_core_count; ++core_index)
    {
        const std::vector<Task>& worker_tasks =
            assignments[static_cast<std::size_t>(core_index - 1)];

        output << "    {";

        for (int column = 0; column < maximum_jobs; ++column)
        {
            if (column > 0)
            {
                output << ", ";
            }

            if (column < static_cast<int>(worker_tasks.size()))
            {
                output << worker_tasks[static_cast<std::size_t>(column)].job_id;
            }
            else
            {
                output << -1;
            }
        }

        output << "}";

        if (core_index < total_core_count - 1)
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n\n";
}


void writeExpectedLoads1(
    std::ofstream& output,
    int total_core_count,
    const std::vector<LoadValue>& worker_loads
)
{
    output
        << "inline constexpr long long expected_load_value["
        << total_core_count
        << "] = {\n";

    output << "    0,\n";

    for (int core_index = 1; core_index < total_core_count; ++core_index)
    {
        output
            << "    "
            << worker_loads[static_cast<std::size_t>(core_index - 1)];

        if (core_index < total_core_count - 1)
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n";
}


void createHeaderFile(
    const std::string& filename,
    int total_core_count,
    const std::vector<Task>& tasks
)
{
    if (total_core_count < 2)
    {
        throw std::invalid_argument(
            "total_core_count must be at least 2."
        );
    }

    const int worker_count = total_core_count - 1;

    if (worker_count > static_cast<int>(tasks.size()))
    {
        throw std::invalid_argument(
            "The calculating worker count exceeds the number of jobs."
        );
    }

    const DistributionMethod method =
        chooseDistributionMethod(worker_count);

    const Assignments assignments =
        distributeTasks(tasks, worker_count, method);

    const DistributionResult result =
        calculateDistributionResult(method, assignments);

    const int maximum_jobs = calculateMaximumJobs(result.assignments);

    std::ofstream output(filename);

    if (!output)
    {
        throw std::runtime_error(
            "Could not create header file: " + filename
        );
    }

    output
        << "#pragma once\n\n"
        << "// Generated automatically by headerCreation.cc\n"
        << "// Total MPI cores: " << total_core_count << "\n"
        << "// Calculating workers: " << worker_count << "\n"
        << "// Distribution method: "
        << distributionMethodName(method) << "\n"
        << "// Maximum expected worker load: "
        << result.maximum_load << "\n\n";

    output
        << "inline constexpr int max_jobs = "
        << maximum_jobs
        << ";\n\n";

    writeCoreJobsMatrix1(
        output,
        total_core_count,
        maximum_jobs,
        result.assignments
    );

    writeExpectedLoads1(
        output,
        total_core_count,
        result.worker_loads
    );

    std::cout
        << "Created "
        << filename
        << " using "
        << distributionMethodName(method)
        << " for "
        << worker_count
        << " workers.\n";
}


// REAL COLLISION-INTEGRAL TASK CREATION

/*
N_bins follows calculate_diagnostics.cc:

    N_bins = N_trap + 5

    estimate_load() is calculated once per physical epsilon bin.
    The neutrino and antineutrino jobs for the same bin receive the
    same estimated load, matching the existing diagnostics workflow.
*/
std::vector<Task> createCollisionIntegralTasks(
    int number_of_trapezoid_bins,
    double epsilon_maximum
)
{
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

    std::cout
        << "Creating epsilon grid and calculating load factors...\n";

    linspace_and_gl* epsilon_grid =
        new linspace_and_gl(
            0.0,
            epsilon_maximum,
            number_of_trapezoid_bins,
            5
        );

    const int number_of_bins =
        number_of_trapezoid_bins + 5;

    std::vector<LoadValue> load_factors(
        static_cast<std::size_t>(number_of_bins),
        0
    );

    collision_integral** integrators =
        new collision_integral*[
            static_cast<std::size_t>(number_of_bins)
        ];

    for (
        int bin = 0;
        bin < number_of_bins;
        ++bin
    )
    {
        integrators[bin] =
            new nu_nu_collision(
                bin,
                epsilon_grid,
                true
            );
    }

    for (
        int bin = 0;
        bin < number_of_bins;
        ++bin
    )
    {
        load_factors[
            static_cast<std::size_t>(bin)
        ] =
            static_cast<LoadValue>(
                integrators[bin]->estimate_load()
            );

        std::cout
            << "\rCalculated load factor "
            << bin + 1
            << " of "
            << number_of_bins
            << std::flush;
    }

    std::cout << "\nLoad-factor calculation finished.\n";

    for (
        int bin = 0;
        bin < number_of_bins;
        ++bin
    )
    {
        delete integrators[bin];
    }

    //delete[] integrators;
    //delete epsilon_grid;
    integrators = nullptr;
    epsilon_grid = nullptr;

    std::vector<Task> tasks;

    tasks.reserve(
        static_cast<std::size_t>(
            2 * number_of_bins
        )
    );

    // Neutrino jobs:
    //
    //     0 through N_bins - 1
    for (
        int bin = 0;
        bin < number_of_bins;
        ++bin
    )
    {
        tasks.push_back({
            bin,
            bin,
            load_factors[
                static_cast<std::size_t>(bin)
            ]
        });
    }

    // Antineutrino jobs:
    //
    //     N_bins through 2*N_bins - 1
    //
    // Therefore:
    //
    //     job_id % N_bins
    //
    // recovers the physical bin, and:
    //
    //     job_id < N_bins
    //
    // is true only for neutrino jobs.
    for (
        int bin = 0;
        bin < number_of_bins;
        ++bin
    )
    {
        tasks.push_back({
            number_of_bins + bin,
            bin,
            load_factors[
                static_cast<std::size_t>(bin)
            ]
        });
    }

    std::cout
        << "Created "
        << tasks.size()
        << " collision-integral jobs from "
        << number_of_bins
        << " physical bins.\n\n";

    return tasks;
}


// OUTPUT-FOLDER HELPERS

void createOutputFolder(
    const std::filesystem::path& output_folder
)
{
    if (output_folder.empty())
    {
        throw std::invalid_argument(
            "The output-folder path cannot be empty."
        );
    }

    std::error_code error;

    std::filesystem::create_directories(
        output_folder,
        error
    );

    if (error)
    {
        throw std::runtime_error(
            "Could not create output folder: "
            + output_folder.string()
            + "\nReason: "
            + error.message()
        );
    }
}


/* Deletes everything inside the selected output folder while
    preserving the output folder itself.

    In main(), comment out:

        clearOutputFolder(output_folder);

whenever the existing generated headers should be kept.
*/
void clearOutputFolder(
    const std::filesystem::path& output_folder
)
{
    createOutputFolder(output_folder);

    std::cout
        << "Deleting the old contents of:\n"
        << "    "
        << output_folder.string()
        << "\n";

    std::error_code error;

    for (
        const std::filesystem::directory_entry& entry
        : std::filesystem::directory_iterator(
            output_folder,
            error
        )
    )
    {
        if (error)
        {
            throw std::runtime_error(
                "Could not read output folder: "
                + output_folder.string()
                + "\nReason: "
                + error.message()
            );
        }

        std::filesystem::remove_all(
            entry.path(),
            error
        );

        if (error)
        {
            throw std::runtime_error(
                "Could not delete: "
                + entry.path().string()
                + "\nReason: "
                + error.message()
            );
        }
    }

    std::cout
        << "The output folder is now empty.\n\n";
}


// HEADER CREATION HELPERS

int calculateMaximumJobs1(
    const Assignments& assignments
)
{
    int maximum_jobs = 0;

    for (
        const std::vector<Task>& worker_tasks
        : assignments
    )
    {
        maximum_jobs = std::max(
            maximum_jobs,
            static_cast<int>(
                worker_tasks.size()
            )
        );
    }

    return maximum_jobs;
}


void writeCoreJobsMatrix(
    std::ofstream& output,
    int total_core_count,
    int maximum_jobs,
    const Assignments& assignments
)
{
    output
        << "inline constexpr int core_jobs["
        << total_core_count
        << "]["
        << maximum_jobs
        << "] = {\n";

    // MPI core 0 is the coordinator.
    output << "    {";

    for (
        int column = 0;
        column < maximum_jobs;
        ++column
    )
    {
        if (column > 0)
        {
            output << ", ";
        }

        output << -1;
    }

    output << "}";

    if (total_core_count > 1)
    {
        output << ",";
    }

    output << "\n";

    // assignments[0] belongs to MPI core 1,
    // assignments[1] belongs to MPI core 2, and so forth.
    for (
        int core_index = 1;
        core_index < total_core_count;
        ++core_index
    )
    {
        const std::vector<Task>& worker_tasks =
            assignments[
                static_cast<std::size_t>(
                    core_index - 1
                )
            ];

        output << "    {";

        for (
            int column = 0;
            column < maximum_jobs;
            ++column
        )
        {
            if (column > 0)
            {
                output << ", ";
            }

            if (
                column
                < static_cast<int>(
                    worker_tasks.size()
                )
            )
            {
                output
                    << worker_tasks[
                        static_cast<std::size_t>(
                            column
                        )
                    ].job_id;
            }
            else
            {
                output << -1;
            }
        }

        output << "}";

        if (core_index < total_core_count - 1)
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n\n";
}


void writeExpectedLoads(
    std::ofstream& output,
    int total_core_count,
    const std::vector<LoadValue>& worker_loads
)
{
    output
        << "inline constexpr long long "
        << "expected_load_value["
        << total_core_count
        << "] = {\n";

    // Core 0 receives no collision-integral work.
    output << "    0";

    if (total_core_count > 1)
    {
        output << ",";
    }

    output << "\n";

    for (
        int core_index = 1;
        core_index < total_core_count;
        ++core_index
    )
    {
        output
            << "    "
            << worker_loads[
                static_cast<std::size_t>(
                    core_index - 1
                )
            ];

        if (core_index < total_core_count - 1)
        {
            output << ",";
        }

        output << "\n";
    }

    output << "};\n";
}


void createHeaderFile(
    const std::filesystem::path& filename,
    int total_core_count,
    int number_of_trapezoid_bins,
    double epsilon_maximum,
    const std::vector<Task>& tasks
)
{
    if (total_core_count < 2)
    {
        throw std::invalid_argument(
            "N_cores must be at least 2."
        );
    }

    const int worker_count =
        total_core_count - 1;

    if (
        worker_count
        > static_cast<int>(
            tasks.size()
        )
    )
    {
        throw std::invalid_argument(
            "The calculating-worker count cannot exceed "
            "the number of collision-integral jobs."
        );
    }

    const DistributionMethod method =
        chooseDistributionMethod(
            worker_count
        );

    std::cout
        << "Selected distribution method: "
        << distributionMethodName(method)
        << "\n"
        << "Calculating workers: "
        << worker_count
        << "\n"
        << "Distributing jobs...\n";

    const Assignments assignments =
        distributeTasks(
            tasks,
            worker_count,
            method
        );

    const DistributionResult result =
        calculateDistributionResult(
            method,
            assignments
        );

    const int maximum_jobs =
        calculateMaximumJobs(
            result.assignments
        );

    std::ofstream output(filename);

    if (!output)
    {
        throw std::runtime_error(
            "Could not create header file: "
            + filename.string()
        );
    }

    const int number_of_bins =
        number_of_trapezoid_bins + 5;

    output
        << "#pragma once\n\n"
        << "// Generated automatically by headerCreation.cc.\n"
        << "// N_cores = "
        << total_core_count
        << "\n"
        << "// N_trap = "
        << number_of_trapezoid_bins
        << "\n"
        << "// eps_max = "
        << std::setprecision(17)
        << epsilon_maximum
        << "\n"
        << "// N_bins = "
        << number_of_bins
        << "\n"
        << "// Calculating workers = "
        << worker_count
        << "\n"
        << "// Distribution method = "
        << distributionMethodName(method)
        << "\n"
        << "// Maximum expected worker load = "
        << result.maximum_load
        << "\n\n";

    output
        << "inline constexpr int max_jobs = "
        << maximum_jobs
        << ";\n\n";

    writeCoreJobsMatrix(
        output,
        total_core_count,
        maximum_jobs,
        result.assignments
    );

    writeExpectedLoads(
        output,
        total_core_count,
        result.worker_loads
    );

    output.close();

    std::cout
        << "Header creation finished.\n"
        << "Maximum jobs assigned to one worker: "
        << maximum_jobs
        << "\n"
        << "Maximum expected worker load: "
        << result.maximum_load
        << "\n";
}


// MAIN


int Main(int argc, char* argv[])
{
    try
    {
        /* Required command:
        
            headerCreation.exe N_cores N_trap eps_max
        
        Example:
        
            headerCreation.exe 128 201 20
        */
       
        if (argc != 4)
        {
            std::cerr
                << "Usage:\n"
                << "    "
                << argv[0]
                << " N_cores N_trap eps_max\n\n"
                << "Example:\n"
                << "    "
                << argv[0]
                << " 128 201 20\n";

            return 1;
        }
        

        const int number_of_cores = std::stoi(argv[1]);

        const int number_of_trapezoid_bins = std::stoi(argv[2]);

        const double epsilon_maximum = std::stod(argv[3]);

        //const int number_of_cores = 128;

        //const int number_of_trapezoid_bins = 201;

        //const int epsilon_maximum = 20;

        if (number_of_cores < 2)
        {
            throw std::invalid_argument(
                "N_cores must be at least 2."
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


        // OUTPUT SETTING

        // Change this one path to choose where the generated
        // header file will be written.
        //
        // Relative example:
        //     "./generated_headers"
        //
        // Windows example:
        //     "C:/Users/user/Desktop/workspace/generated_headers"
        const std::filesystem::path output_folder =
            "../../WorkerHeaderCreation/headers";
            //"./generated_headers";
            //creating file to use as the outputfolder


        std::cout
            << "Starting headerCreation.cc\n"
            << "N_cores: "
            << number_of_cores
            << "\n"
            << "N_trap: "
            << number_of_trapezoid_bins
            << "\n"
            << "eps_max: "
            << epsilon_maximum
            << "\n"
            << "Output folder:\n"
            << "    "
            << output_folder.string()
            << "\n\n";


        // This removes everything currently inside output_folder.
        //
        // COMMENT OUT THIS ONE LINE whenever the old files should
        // remain in the folder.
        clearOutputFolder(output_folder);


        // This is still necessary when clearOutputFolder() is
        // commented out.
        createOutputFolder(output_folder);


        const std::vector<Task> tasks =
            createCollisionIntegralTasks(
                number_of_trapezoid_bins,
                epsilon_maximum
            );


        const std::filesystem::path output_filename =
            output_folder
            / (
                "collision_distribution_"
                + std::to_string(number_of_cores)
                + "_cores_"
                + std::to_string(number_of_trapezoid_bins)
                + "_trap.hh"
            );


        createHeaderFile(
            output_filename,
            number_of_cores,
            number_of_trapezoid_bins,
            epsilon_maximum,
            tasks
        );


        std::cout
            << "\nFinished successfully.\n"
            << "Created:\n"
            << "    "
            << std::filesystem::absolute(
                output_filename
            ).string()
            << "\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "\nheaderCreation.cc stopped because of an error:\n"
            << error.what()
            << "\n";

        return 1;
    }
}



/*
task:
    create a header file with the following:
        int max_jobs = ##;
        int core_jobs[N_cores][##] = {…};
        Int expected_load_value[N_cores] = {…};
    
    intputs: 
        number of cores, 
        number of trapezoid bins, 
        and eps_max
*/

/*
    currently have the distribution methods.
    now i need to go ahead and form a way to 
*/

/*
we need to output a 1 file, that one file will will have the 4
    epsilon values
    the trapazoid bins
    the load factors
    all from an output from a file. which is calculate_diagnostics.cc 

*/

/*
Here’s what we want: choose a method to distribute the collision integrals, and then we need a script that does the following:

Inputs: (int) N_cores, (int) N_trap, (double) eps_max

(1) given (int) N_trap and (double) eps_max, we want to create c++ code that determines all of the load factors for each bin. Test with N_trap = 201 and eps_max = 20 which gives us the load factors from the data file you already have.
See run/calculate_diagnostics.cc for using the linspace_and_gl and collision_integral classes

(2) find an optimal strategy to use N_cores. Print a header file that has:
int max_jobs = ##;
int core_jobs[N_cores][##] = {…};
Int expected_load_value[N_cores] = {…};

The matrix core_jobs should have the following properties:
- The number of columns should equal the maximum jobs for a single core
- Every row of core_jobs—core_jobs[#][…]—should include all the collision integrals for that specific worker to calculate. The numbering scheme should be an integer between 0 and 2*N_bins, with integers greater than N_bins recognized as antineutrino integrals (so that "% N_bins" can be used to identify the specific bin, and "< N_bins" can be used to determine if it is neutrino (true) or antineutrino (false)
- Whenever a core has less than the maximum jobs, use -1 in the remaining values
- The first row, core_jobs[0][…], should all be -1 because worker 0 does not calculate collision integrals, no matter what.

The expected_load_value array is there to double check.

Overall, your script should accept three command line inputs: number of cores, number of trapezoid bins, and eps_max
Your script should create the header file. If any intermediate files are created, they should be deleted.

(1) needs to be written in c++
(2) can either interface with (1) via text file; or I guess you could c++ it too and put it with (1)

*/