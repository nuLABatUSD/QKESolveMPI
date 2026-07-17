import csv
import heapq
import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np


filename = "R_values-201-20-5.csv"
#output_folder = "./solvingProblems/workerFolder"


INPUT_FILE = "./solvingProblems/" + filename
OUTPUT_FOLDER = "./solvingProblems/workload_comparison"

NUMBER_OF_WORKERS = 127

NUMBER_OF_EPSILON_BINS = 205

# Four rows are expected:
#   row 0: R
#   row 1: epsilon
#   row 2: quadrature weights
#   row 3: load factors


@dataclass
class Task:
    """
    One collision-integral calculation.

    task_number:
        Global task number from 1 to 410.

    particle_type:
        "neutrino" or "antineutrino".

    epsilon_bin:
        Epsilon-bin number from 1 to 205.

    load:
        Estimated calculation cost.
    """

    task_number: int
    particle_type: str
    epsilon_bin: int
    load: float

def read_four_row_file(filename):
    """
    Reads either a comma-separated or whitespace-separated file.

    Returns:
        R_values
        epsilon_values
        weight_values
        load_factors
    """

    rows = []

    with open(filename, "r", newline="") as file:
        for raw_line in file:
            raw_line = raw_line.strip()

            if not raw_line:
                continue

            if "," in raw_line:
                values = next(csv.reader([raw_line]))
            else:
                values = raw_line.split()

            cleaned = [
                float(value.strip())
                for value in values
                if value.strip() != ""
            ]

            rows.append(cleaned)

    if len(rows) != 4:
        raise ValueError(
            f"Expected exactly 4 nonempty rows, but found {len(rows)}."
        )

    row_lengths = [len(row) for row in rows]

    if any(length != 206 for length in row_lengths):
        raise ValueError(
            "Expected 206 columns in every row.\n"
            f"Found row lengths: {row_lengths}"
        )

    return (
        np.asarray(rows[0], dtype=float),
        np.asarray(rows[1], dtype=float),
        np.asarray(rows[2], dtype=float),
        np.asarray(rows[3], dtype=float),
    )

def create_tasks(load_factors):
    """
    Creates the 410 tasks:

        Tasks 1-205:
            neutrino bins 1-205

        Tasks 206-410:
            antineutrino bins 1-205

    Epsilon bin 0 is ignored.
    """

    tasks = []

    task_number = 1

    for particle_type in ("neutrino", "antineutrino"):
        for epsilon_bin in range(1, NUMBER_OF_EPSILON_BINS + 1):
            tasks.append(
                Task(
                    task_number=task_number,
                    particle_type=particle_type,
                    epsilon_bin=epsilon_bin,
                    load=float(load_factors[epsilon_bin]),
                )
            )

            task_number += 1

    return tasks

def create_empty_assignments(worker_count):
    return [[] for _ in range(worker_count)]


def calculate_worker_loads(assignments):
    """
    Adds the loads of all tasks assigned to each worker.
    """

    return np.asarray(
        [
            sum(task.load for task in worker_tasks)
            for worker_tasks in assignments
        ],
        dtype=float,
    )


def calculate_statistics(assignments):
    worker_loads = calculate_worker_loads(assignments)

    average_load = float(np.mean(worker_loads))
    maximum_load = float(np.max(worker_loads))
    minimum_load = float(np.min(worker_loads))
    standard_deviation = float(np.std(worker_loads))
    load_range = maximum_load - minimum_load

    if maximum_load > 0:
        efficiency = average_load / maximum_load
    else:
        efficiency = 0.0

    return {
        "loads": worker_loads,
        "average": average_load,
        "maximum": maximum_load,
        "minimum": minimum_load,
        "range": load_range,
        "standard_deviation": standard_deviation,
        "efficiency": efficiency,
        "largest_worker": int(np.argmax(worker_loads)) + 1,
        "smallest_worker": int(np.argmin(worker_loads)) + 1,
    }


# METHOD 1: WALL


def distribute_wall(tasks, worker_count):
    """
    ordinary role assignment.

    example with 127 workers:

        tasks 1-127   -> workers 1-127
        tasks 128-254 -> workers 1-127
        tasks 255-381 -> workers 1-127
        tasks 382-410 -> workers 1-29
    """

    assignments = create_empty_assignments(worker_count)

    for task_index, task in enumerate(tasks):
        worker_index = task_index % worker_count
        assignments[worker_index].append(task)

    return assignments


# METHOD 2: SNAKE


def distribute_snake(tasks, worker_count):
    """
    alternates assignment direction every complete or partial row.

    first row:
        worker 1 -> worker N

    second row:
        worker N -> worker 1

    third row:
        worker 1 -> worker N

    and so on.
    """

    assignments = create_empty_assignments(worker_count)

    for task_index, task in enumerate(tasks):
        row_number = task_index // worker_count
        position_in_row = task_index % worker_count

        if row_number % 2 == 0:
            worker_index = position_in_row
        else:
            worker_index = worker_count - 1 - position_in_row

        assignments[worker_index].append(task)

    return assignments


# METHOD 3: MOUNTAIN / CONTIGUOUS MINIMAX PARTITION


def number_of_groups_needed(tasks, maximum_group_load):
    """
    determines how many contiguous groups are required if no group
    may exceed maximum_group_load.

    this assumes all task loads are nonnegative.
    """

    group_count = 1
    current_load = 0.0

    for task in tasks:
        if task.load > maximum_group_load:
            return len(tasks) + 1

        if current_load + task.load <= maximum_group_load:
            current_load += task.load
        else:
            group_count += 1
            current_load = task.load

    return group_count


def find_contiguous_load_limit(tasks, worker_count):
    """
    binary searches for the smallest maximum load that allows the
    ordered task sequence to be split into no more than worker_count
    contiguous groups.

    this is a more systematic version of the proposed mountain method.
    """

    lower_bound = max(task.load for task in tasks)
    upper_bound = sum(task.load for task in tasks)

    # integer load factors are expected, but float arithmetic is allowed.
    # sixty iterations is sufficient for high numerical precision.
    for _ in range(60):
        middle = (lower_bound + upper_bound) / 2.0

        required_groups = number_of_groups_needed(tasks, middle)

        if required_groups <= worker_count:
            upper_bound = middle
        else:
            lower_bound = middle

    return upper_bound


def build_exact_contiguous_groups(tasks, worker_count, load_limit):
    """
    divides the ordered task list into exactly worker_count nonempty,
    contiguous groups.

    different workers may receive different numbers of tasks.
    The goal is to keep every worker's total load at or below load_limit
    whenever possible.
    """

    if worker_count <= 0:
        raise ValueError("worker_count must be greater than zero.")

    if worker_count > len(tasks):
        raise ValueError(
            "Cannot create more nonempty worker groups than tasks.\n"
            f"Workers: {worker_count}\n"
            f"Tasks: {len(tasks)}"
        )

    assignments = []
    current_group = []
    current_load = 0.0

    for task_index, task in enumerate(tasks):
        tasks_remaining_including_current = len(tasks) - task_index

        # If we finalize current_group now, this is how many workers
        # would still need to receive tasks.
        workers_remaining_after_current_group = (
            worker_count - len(assignments) - 1
        )

        # We must close the current group before this task when the
        # number of remaining tasks exactly equals the number of
        # remaining workers. Each remaining worker must then receive
        # at least one task.
        must_start_new_group = (
            len(current_group) > 0
            and tasks_remaining_including_current
            == workers_remaining_after_current_group
        )

        exceeds_limit = (
            len(current_group) > 0
            and current_load + task.load > load_limit
        )

        # Do not create more than worker_count - 1 completed groups.
        can_start_new_group = (
            len(assignments) < worker_count - 1
        )

        if can_start_new_group and (
            must_start_new_group or exceeds_limit
        ):
            assignments.append(current_group)
            current_group = []
            current_load = 0.0

        current_group.append(task)
        current_load += task.load

    if current_group:
        assignments.append(current_group)

    if len(assignments) != worker_count:
        raise RuntimeError(
            "Mountain partition did not produce exactly "
            f"{worker_count} groups. Produced {len(assignments)}.\n"
            f"Load limit: {load_limit}"
        )

    return assignments

def distribute_mountain(tasks, worker_count):
    """
    keeps the 410 tasks in their original order and finds a contiguous
    partition that minimizes the largest worker load.

    workers near the beginning may receive several inexpensive tasks,
    while workers near the end may receive fewer expensive tasks.
    """

    load_limit = find_contiguous_load_limit(tasks, worker_count)

    return build_exact_contiguous_groups(
        tasks,
        worker_count,
        load_limit,
    )


# METHOD 4: LARGEST-PROCESSING-TIME GREEDY


def distribute_largest_first(tasks, worker_count):
    """
    sorts tasks from largest load to smallest load.

    each task is assigned to the worker with the smallest current total.
    A min-heap makes this efficient.
    """

    assignments = create_empty_assignments(worker_count)

    # heap entries:
    # current_worker_load, worker_index
    worker_heap = [
        (0.0, worker_index)
        for worker_index in range(worker_count)
    ]

    heapq.heapify(worker_heap)

    sorted_tasks = sorted(
        tasks,
        key=lambda task: task.load,
        reverse=True,
    )

    for task in sorted_tasks:
        current_load, worker_index = heapq.heappop(worker_heap)

        assignments[worker_index].append(task)

        new_load = current_load + task.load

        heapq.heappush(
            worker_heap,
            (new_load, worker_index),
        )

    return assignments


# OPTIONAL METHOD 5: LARGEST-SMALLEST PAIRING

def distribute_high_low(tasks, worker_count):
    """
    creates a task ordering that alternates between the largest and
    smallest remaining task, then applies snake assignment.

    this is included as another inexpensive heuristic.
    """

    sorted_tasks = sorted(tasks, key=lambda task: task.load)

    reordered_tasks = []

    left = 0
    right = len(sorted_tasks) - 1

    while left <= right:
        reordered_tasks.append(sorted_tasks[right])
        right -= 1

        if left <= right:
            reordered_tasks.append(sorted_tasks[left])
            left += 1

    return distribute_snake(reordered_tasks, worker_count)


# END OF ALL METHOD TYPES ====================================

def save_assignments_csv(method_name, assignments, output_folder):
    output_path = os.path.join(
        output_folder,
        f"{method_name}_assignments.csv",
    )

    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "worker",
            "number_of_tasks",
            "total_load",
            "task_numbers",
            "particle_types",
            "epsilon_bins",
            "individual_loads",
        ])

        for worker_index, worker_tasks in enumerate(assignments):
            writer.writerow([
                worker_index + 1,
                len(worker_tasks),
                sum(task.load for task in worker_tasks),
                " ".join(
                    str(task.task_number)
                    for task in worker_tasks
                ),
                " ".join(
                    task.particle_type
                    for task in worker_tasks
                ),
                " ".join(
                    str(task.epsilon_bin)
                    for task in worker_tasks
                ),
                " ".join(
                    str(task.load)
                    for task in worker_tasks
                ),
            ])


def save_summary_csv(results, output_folder, lower_bound):
    output_path = os.path.join(
        output_folder,
        "method_summary.csv",
    )

    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "method",
            "maximum_load",
            "minimum_load",
            "average_load",
            "load_range",
            "standard_deviation",
            "balance_efficiency",
            "largest_worker",
            "smallest_worker",
            "theoretical_lower_bound",
            "maximum_above_lower_bound",
        ])

        for method_name, result in results.items():
            writer.writerow([
                method_name,
                result["maximum"],
                result["minimum"],
                result["average"],
                result["range"],
                result["standard_deviation"],
                result["efficiency"],
                result["largest_worker"],
                result["smallest_worker"],
                lower_bound,
                result["maximum"] - lower_bound,
            ])


# PLOTS


def plot_worker_loads(results, output_folder):
    """
    shows the methods in their actual worker-number order.
    """

    plt.figure(figsize=(15, 8))

    for method_name, result in results.items():
        worker_numbers = np.arange(1, len(result["loads"]) + 1)

        plt.plot(
            worker_numbers,
            result["loads"],
            marker="o",
            markersize=2.5,
            linewidth=1,
            label=(
                f"{method_name}: "
                f"maximum = {result['maximum']:.0f}"
            ),
        )

    average_load = next(iter(results.values()))["average"]

    plt.axhline(
        average_load,
        linestyle="--",
        linewidth=1,
        label=f"Average load = {average_load:.2f}",
    )

    plt.xlabel("Worker number")
    plt.ylabel("Total load factor")
    plt.title(
        "Collision-Integral Workload Distribution\n"
        f"{NUMBER_OF_WORKERS} workers, 410 collision integrals"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "all_methods_worker_order.png",
    )

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_sorted_worker_loads(results, output_folder):
    """
    sorts each method's worker loads from smallest to largest.

    this removes arbitrary worker numbering and makes the degree of
    imbalance easier to compare.
    """

    plt.figure(figsize=(15, 8))

    for method_name, result in results.items():
        sorted_loads = np.sort(result["loads"])
        rank = np.arange(1, len(sorted_loads) + 1)

        plt.plot(
            rank,
            sorted_loads,
            marker="o",
            markersize=2.5,
            linewidth=1,
            label=(
                f"{method_name}: "
                f"range = {result['range']:.0f}"
            ),
        )

    average_load = next(iter(results.values()))["average"]

    plt.axhline(
        average_load,
        linestyle="--",
        linewidth=1,
        label=f"Average load = {average_load:.2f}",
    )

    plt.xlabel("Worker load rank")
    plt.ylabel("Total load factor")
    plt.title(
        "Sorted Worker Loads\n"
        "Flatter curves indicate better load balance"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "all_methods_sorted_loads.png",
    )

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_maximum_loads(results, output_folder, lower_bound):
    method_names = list(results.keys())

    maximum_loads = [
        results[method_name]["maximum"]
        for method_name in method_names
    ]

    plt.figure(figsize=(12, 7))

    bars = plt.bar(
        method_names,
        maximum_loads,
    )

    plt.axhline(
        lower_bound,
        linestyle="--",
        linewidth=1,
        label=f"Theoretical lower bound = {lower_bound:.2f}",
    )

    for bar, maximum_load in zip(bars, maximum_loads):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{maximum_load:.0f}",
            ha="center",
            va="bottom",
        )

    plt.xlabel("Distribution method")
    plt.ylabel("Largest worker load")
    plt.title(
        "Largest Worker Load by Distribution Method\n"
        "Lower is better"
    )

    plt.grid(True, axis="y", alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "maximum_load_comparison.png",
    )

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_task_counts(results, assignments_by_method, output_folder):
    """
    shows how many collision integrals each worker receives.
    """

    plt.figure(figsize=(15, 8))

    for method_name, assignments in assignments_by_method.items():
        task_counts = [
            len(worker_tasks)
            for worker_tasks in assignments
        ]

        workers = np.arange(1, len(task_counts) + 1)

        plt.plot(
            workers,
            task_counts,
            marker="o",
            markersize=2.5,
            linewidth=1,
            label=method_name,
        )

    plt.xlabel("Worker number")
    plt.ylabel("Number of collision integrals")
    plt.title("Number of Tasks Assigned to Each Worker")

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "worker_task_counts.png",
    )

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

# PRINT RESULTS

def print_results(results, lower_bound):
    """
    prints the comparison of all distribution methods for one
    particular worker count.
    """

    print()
    print("=" * 100)
    print("WORKLOAD DISTRIBUTION RESULTS")
    print("=" * 100)

    heading = (
        f"{'Method':<25}"
        f"{'Maximum':>14}"
        f"{'Minimum':>14}"
        f"{'Average':>14}"
        f"{'Range':>14}"
        f"{'Std. dev.':>14}"
        f"{'Efficiency':>14}"
    )

    print(heading)
    print("-" * len(heading))

    sorted_results = sorted(
        results.items(),
        key=lambda item: item[1]["maximum"],
    )

    for method_name, result in sorted_results:
        print(
            f"{method_name:<25}"
            f"{result['maximum']:>14.2f}"
            f"{result['minimum']:>14.2f}"
            f"{result['average']:>14.2f}"
            f"{result['range']:>14.2f}"
            f"{result['standard_deviation']:>14.2f}"
            f"{100.0 * result['efficiency']:>13.2f}%"
        )

    best_method, best_result = sorted_results[0]

    print()
    print(f"Theoretical lower bound: {lower_bound:.2f}")
    print(f"Best tested method:       {best_method}")
    print(f"Best maximum load:        {best_result['maximum']:.2f}")

    print(
        "Distance above lower bound: "
        f"{best_result['maximum'] - lower_bound:.2f}"
    )

    print()
    print(
        f"Slowest worker for {best_method}: "
        f"worker {best_result['largest_worker']}"
    )

    print(
        f"Fastest worker for {best_method}: "
        f"worker {best_result['smallest_worker']}"
    )
    
def print_reduced_worker_comparison(
    best_by_worker_count,
    reference_worker_count=127,
    minimum_slowdown_percent=1,
    maximum_slowdown_percent=30,
):
    """
    prints the smallest worker count that stays within each allowed
    slowdown percentage relative to the reference worker count.
    """

    if reference_worker_count not in best_by_worker_count:
        raise ValueError(
            f"Reference worker count {reference_worker_count} "
            "was not included in the benchmark."
        )

    reference_result = best_by_worker_count[
        reference_worker_count
    ]

    reference_maximum = reference_result["maximum"]
    reference_method = reference_result["method"]

    print()
    print("=" * 100)
    print("REDUCED WORKER COUNT COMPARISON")
    print("=" * 100)

    print(
        f"Reference workers:      "
        f"{reference_worker_count}"
    )

    print(
        f"Reference best method:  "
        f"{reference_method}"
    )

    print(
        f"Reference maximum load: "
        f"{reference_maximum:.2f}"
    )

    print("-" * 100)

    for allowed_slowdown in range(
        minimum_slowdown_percent,
        maximum_slowdown_percent + 1,
    ):
        reduced_result = (
            find_smallest_worker_count_near_reference(
                best_by_worker_count,
                reference_worker_count=reference_worker_count,
                allowed_slowdown_percent=float(
                    allowed_slowdown
                ),
            )
        )

        if reduced_result is None:
            print(
                f"Within {allowed_slowdown:2d}%: "
                "No valid worker count found."
            )
            continue

        reduced_worker_count = reduced_result["worker_count"]

        workers_removed = (
            reference_worker_count
            - reduced_worker_count
        )

        worker_reduction_percent = (
            workers_removed
            / reference_worker_count
            * 100.0
        )

        actual_slowdown_percent = (
            reduced_result["maximum"]
            / reference_maximum
            - 1.0
        ) * 100.0

        print(
            f"Within {allowed_slowdown:2d}%: "
            f"{reduced_worker_count:3d} workers, "
            f"remove {workers_removed:3d} "
            f"({worker_reduction_percent:6.2f}%), "
            f"using {reduced_result['method']:<18} "
            f"maximum load "
            f"{reduced_result['maximum']:12.2f}, "
            f"actual slowdown "
            f"{actual_slowdown_percent:6.2f}%"
        )    
        
def save_reduced_worker_comparison(
    best_by_worker_count,
    reference_worker_count,
    minimum_slowdown_percent,
    maximum_slowdown_percent,
    output_folder,
):
    output_path = os.path.join(
        output_folder,
        "reduced_worker_comparison.csv",
    )

    reference_maximum = (
        best_by_worker_count[reference_worker_count]["maximum"]
    )

    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "allowed_slowdown_percent",
            "worker_count",
            "workers_removed",
            "worker_reduction_percent",
            "best_method",
            "maximum_load",
            "reference_maximum_load",
            "actual_slowdown_percent",
            "allowed_maximum_load",
        ])

        for allowed_slowdown in range(
            minimum_slowdown_percent,
            maximum_slowdown_percent + 1,
        ):
            result = find_smallest_worker_count_near_reference(
                best_by_worker_count,
                reference_worker_count=reference_worker_count,
                allowed_slowdown_percent=float(
                    allowed_slowdown
                ),
            )

            if result is None:
                continue

            workers_removed = (
                reference_worker_count
                - result["worker_count"]
            )

            worker_reduction_percent = (
                workers_removed
                / reference_worker_count
                * 100.0
            )

            actual_slowdown_percent = (
                result["maximum"]
                / reference_maximum
                - 1.0
            ) * 100.0

            writer.writerow([
                allowed_slowdown,
                result["worker_count"],
                workers_removed,
                worker_reduction_percent,
                result["method"],
                result["maximum"],
                reference_maximum,
                actual_slowdown_percent,
                result["allowed_maximum"],
            ])

def run_all_distribution_methods(tasks, worker_count):
    """
    runs every distribution method for one specified worker count.

    returns:
        assignments_by_method
        results
        theoretical_lower_bound
    """

    assignments_by_method = {
        "Wall": distribute_wall(
            tasks,
            worker_count,
        ),

        "Snake": distribute_snake(
            tasks,
            worker_count,
        ),

        "Mountain": distribute_mountain(
            tasks,
            worker_count,
        ),

        "Largest-first": distribute_largest_first(
            tasks,
            worker_count,
        ),

        "High-low snake": distribute_high_low(
            tasks,
            worker_count,
        ),
    }

    results = {
        method_name: calculate_statistics(assignments)
        for method_name, assignments
        in assignments_by_method.items()
    }

    total_load = sum(task.load for task in tasks)
    largest_single_task = max(task.load for task in tasks)

    average_load = total_load / worker_count

    theoretical_lower_bound = max(
        average_load,
        largest_single_task,
    )

    return (
        assignments_by_method,
        results,
        theoretical_lower_bound,
    )
    
def benchmark_worker_counts(tasks, minimum_workers, maximum_workers):
    """
    tests every distribution method for every worker count in the
    requested range.

    example:
        minimum_workers = 1
        maximum_workers = 127

    returns a nested dictionary:

        benchmark_results[worker_count][method_name]

    each method result contains its maximum, minimum, average,
    standard deviation, efficiency, and related statistics.
    """

    if minimum_workers < 1:
        raise ValueError("minimum_workers must be at least 1.")

    if maximum_workers > len(tasks):
        raise ValueError(
            "maximum_workers cannot exceed the number of tasks."
        )

    if minimum_workers > maximum_workers:
        raise ValueError(
            "minimum_workers cannot be greater than maximum_workers."
        )

    benchmark_results = {}

    for worker_count in range(
        minimum_workers,
        maximum_workers + 1,
    ):
        print(
            f"Testing worker count "
            f"{worker_count}/{maximum_workers}..."
        )

        (
            assignments_by_method,
            method_results,
            theoretical_lower_bound,
        ) = run_all_distribution_methods(
            tasks,
            worker_count,
        )

        benchmark_results[worker_count] = {
            "methods": method_results,
            "lower_bound": theoretical_lower_bound,
        }

    return benchmark_results

def save_worker_count_benchmark_csv(
    benchmark_results,
    output_folder,
):
    """
    saves one row for every combination of worker count and method.
    """

    output_path = os.path.join(
        output_folder,
        "all_worker_count_results.csv",
    )

    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "worker_count",
            "method",
            "maximum_load",
            "minimum_load",
            "average_load",
            "load_range",
            "standard_deviation",
            "balance_efficiency",
            "theoretical_lower_bound",
            "maximum_above_lower_bound",
            "relative_gap_percent",
            "largest_worker",
            "smallest_worker",
        ])

        for worker_count, worker_data in benchmark_results.items():
            lower_bound = worker_data["lower_bound"]

            for method_name, result in worker_data["methods"].items():
                difference = result["maximum"] - lower_bound

                if lower_bound > 0:
                    relative_gap_percent = (
                        difference / lower_bound
                    ) * 100.0
                else:
                    relative_gap_percent = 0.0

                writer.writerow([
                    worker_count,
                    method_name,
                    result["maximum"],
                    result["minimum"],
                    result["average"],
                    result["range"],
                    result["standard_deviation"],
                    result["efficiency"],
                    lower_bound,
                    difference,
                    relative_gap_percent,
                    result["largest_worker"],
                    result["smallest_worker"],
                ])

def find_best_method_for_each_worker_count(
    benchmark_results,
):
    """
    finds the method with the smallest maximum worker load for every
    worker count.

    returns:
        best_by_worker_count
    """

    best_by_worker_count = {}

    for worker_count, worker_data in benchmark_results.items():
        methods = worker_data["methods"]

        best_method_name, best_result = min(
            methods.items(),
            key=lambda item: item[1]["maximum"],
        )

        best_by_worker_count[worker_count] = {
            "method": best_method_name,
            "maximum": best_result["maximum"],
            "average": best_result["average"],
            "efficiency": best_result["efficiency"],
            "lower_bound": worker_data["lower_bound"],
        }

    return best_by_worker_count


def save_best_method_csv(
    best_by_worker_count,
    output_folder,
    reference_worker_count=127,
):
    """
    saves the best distribution method for every worker count.

    runtime is assumed proportional to the largest worker load:

        runtime_percent =
            current_maximum_load / reference_maximum_load * 100

    cost is assumed proportional to the total number of CPUs,
    including one coordinator:

        worker_cost_percent =
            (current_workers + 1)
            / (reference_workers + 1)
            * 100
    """

    if reference_worker_count not in best_by_worker_count:
        raise ValueError(
            f"Reference worker count {reference_worker_count} "
            "was not included in best_by_worker_count."
        )

    reference_result = best_by_worker_count[
        reference_worker_count
    ]

    reference_maximum_load = reference_result["maximum"]
    reference_total_cores = reference_worker_count + 1

    output_path = os.path.join(
        output_folder,
        "best_method_by_worker_count.csv",
    )

    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "worker_count",
            "total_cores_including_coordinator",
            "best_method",
            "best_maximum_load",
            "average_load",
            "balance_efficiency",
            "theoretical_lower_bound",
            "maximum_above_lower_bound",

            # New columns
            "estimated_runtime_percent",
            "estimated_slowdown_percent",
            "worker_cost_percent",
            "estimated_cost_savings_percent",
        ])

        for worker_count, result in best_by_worker_count.items():
            total_cores = worker_count + 1

            estimated_runtime_percent = (
                result["maximum"]
                / reference_maximum_load
                * 100.0
            )

            estimated_slowdown_percent = (
                estimated_runtime_percent - 100.0
            )

            worker_cost_percent = (
                total_cores
                / reference_total_cores
                * 100.0
            )

            estimated_cost_savings_percent = (
                100.0 - worker_cost_percent
            )

            writer.writerow([
                worker_count,
                total_cores,
                result["method"],
                result["maximum"],
                result["average"],
                result["efficiency"],
                result["lower_bound"],
                result["maximum"] - result["lower_bound"],

                estimated_runtime_percent,
                estimated_slowdown_percent,
                worker_cost_percent,
                estimated_cost_savings_percent,
            ])
            
def count_method_wins(
    benchmark_results,
    tolerance=1e-9,
):
    """
    counts how many worker counts each method wins.

    tied methods are each given one tied win.
    """

    method_win_counts = {}
    ties_by_worker_count = {}

    for worker_count, worker_data in benchmark_results.items():
        methods = worker_data["methods"]

        best_maximum = min(
            result["maximum"]
            for result in methods.values()
        )

        tolerance_value = tolerance * max(
            1.0,
            abs(best_maximum),
        )

        winning_methods = [
            method_name
            for method_name, result in methods.items()
            if abs(
                result["maximum"] - best_maximum
            ) <= tolerance_value
        ]

        ties_by_worker_count[worker_count] = winning_methods

        for method_name in winning_methods:
            method_win_counts[method_name] = (
                method_win_counts.get(method_name, 0) + 1
            )

    return method_win_counts, ties_by_worker_count
            
def save_method_win_counts(
    method_win_counts,
    output_folder,
):
    output_path = os.path.join(
        output_folder,
        "method_win_counts.csv",
    )

    sorted_counts = sorted(
        method_win_counts.items(),
        key=lambda item: item[1],
        reverse=True,
    )

    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "method",
            "number_of_worker_counts_won_or_tied",
        ])

        for method_name, win_count in sorted_counts:
            writer.writerow([
                method_name,
                win_count,
            ])

def plot_maximum_load_vs_worker_count(
    benchmark_results,
    output_folder,
):
    worker_counts = np.asarray(
        list(benchmark_results.keys()),
        dtype=int,
    )

    first_worker_data = next(
        iter(benchmark_results.values())
    )

    method_names = list(
        first_worker_data["methods"].keys()
    )

    plt.figure(figsize=(15, 8))

    for method_name in method_names:
        maximum_loads = np.asarray(
            [
                benchmark_results[worker_count]
                ["methods"][method_name]["maximum"]
                for worker_count in worker_counts
            ],
            dtype=float,
        )

        plt.plot(
            worker_counts,
            maximum_loads,
            linewidth=1.5,
            label=method_name,
        )

    lower_bounds = np.asarray(
        [
            benchmark_results[worker_count]["lower_bound"]
            for worker_count in worker_counts
        ],
        dtype=float,
    )

    plt.plot(
        worker_counts,
        lower_bounds,
        linestyle="--",
        linewidth=1.5,
        label="Theoretical lower bound",
    )

    plt.xlabel("Number of workers")
    plt.ylabel("Largest worker load")
    plt.title(
        "Largest Worker Load Versus Number of Workers\n"
        "Lower values indicate shorter theoretical runtime"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "maximum_load_vs_worker_count.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()  
    
def plot_relative_gap_vs_worker_count(
    benchmark_results,
    output_folder,
):
    worker_counts = np.asarray(
        list(benchmark_results.keys()),
        dtype=int,
    )

    first_worker_data = next(
        iter(benchmark_results.values())
    )

    method_names = list(
        first_worker_data["methods"].keys()
    )

    plt.figure(figsize=(15, 8))

    for method_name in method_names:
        relative_gaps = []

        for worker_count in worker_counts:
            worker_data = benchmark_results[worker_count]
            lower_bound = worker_data["lower_bound"]

            maximum_load = (
                worker_data["methods"]
                [method_name]["maximum"]
            )

            if lower_bound > 0:
                relative_gap = (
                    (maximum_load - lower_bound)
                    / lower_bound
                    * 100.0
                )
            else:
                relative_gap = 0.0

            relative_gaps.append(relative_gap)

        plt.plot(
            worker_counts,
            relative_gaps,
            linewidth=1.5,
            label=method_name,
        )

    plt.axhline(
        0.0,
        linestyle="--",
        linewidth=1,
        label="Perfect theoretical balance",
    )

    plt.xlabel("Number of workers")
    plt.ylabel("Maximum above lower bound (%)")
    plt.title(
        "Distribution Quality Versus Number of Workers\n"
        "Lower percentages indicate better load balancing"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "relative_gap_vs_worker_count.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

def plot_efficiency_vs_worker_count(
    benchmark_results,
    output_folder,
):
    worker_counts = np.asarray(
        list(benchmark_results.keys()),
        dtype=int,
    )

    first_worker_data = next(
        iter(benchmark_results.values())
    )

    method_names = list(
        first_worker_data["methods"].keys()
    )

    plt.figure(figsize=(15, 8))

    for method_name in method_names:
        efficiencies = np.asarray(
            [
                benchmark_results[worker_count]
                ["methods"][method_name]["efficiency"]
                * 100.0
                for worker_count in worker_counts
            ],
            dtype=float,
        )

        plt.plot(
            worker_counts,
            efficiencies,
            linewidth=1.5,
            label=method_name,
        )

    plt.axhline(
        100.0,
        linestyle="--",
        linewidth=1,
        label="Perfect balance",
    )

    plt.xlabel("Number of workers")
    plt.ylabel("Balance efficiency (%)")
    plt.title(
        "Load-Balance Efficiency Versus Number of Workers\n"
        "Higher percentages indicate better balance"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "efficiency_vs_worker_count.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()    
    
def plot_best_tested_load_vs_worker_count(
    best_by_worker_count,
    output_folder,
):
    worker_counts = np.asarray(
        list(best_by_worker_count.keys()),
        dtype=int,
    )

    best_maximums = np.asarray(
        [
            best_by_worker_count[worker_count]["maximum"]
            for worker_count in worker_counts
        ],
        dtype=float,
    )

    lower_bounds = np.asarray(
        [
            best_by_worker_count[worker_count]["lower_bound"]
            for worker_count in worker_counts
        ],
        dtype=float,
    )

    plt.figure(figsize=(15, 8))

    plt.plot(
        worker_counts,
        best_maximums,
        linewidth=1.5,
        label="Best tested distribution",
    )

    plt.plot(
        worker_counts,
        lower_bounds,
        linestyle="--",
        linewidth=1.5,
        label="Theoretical lower bound",
    )

    plt.xlabel("Number of workers")
    plt.ylabel("Largest worker load")
    plt.title(
        "Best Tested Maximum Load Versus Worker Count"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "best_tested_load_vs_worker_count.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()                              

def find_smallest_worker_count_near_reference(
    best_by_worker_count,
    reference_worker_count,
    allowed_slowdown_percent,
):
    """
    finds the smallest worker count whose best tested maximum load
    is no more than the permitted percentage above the reference.
    this asks:
        what is the smallest worker count that is no more than 5%
        slower than the best result at 127 workers?
    """

    if reference_worker_count not in best_by_worker_count:
        raise ValueError(
            "Reference worker count was not benchmarked."
        )

    reference_result = (
        best_by_worker_count[reference_worker_count]
    )

    reference_maximum = reference_result["maximum"]

    allowed_maximum = reference_maximum * (
        1.0 + allowed_slowdown_percent / 100.0
    )

    matching_counts = []

    for worker_count, result in best_by_worker_count.items():
        if result["maximum"] <= allowed_maximum:
            matching_counts.append(worker_count)

    if not matching_counts:
        return None

    smallest_count = min(matching_counts)

    return {
        "worker_count": smallest_count,
        "method": best_by_worker_count[smallest_count]["method"],
        "maximum": best_by_worker_count[smallest_count]["maximum"],
        "reference_maximum": reference_maximum,
        "allowed_maximum": allowed_maximum,
    }

def plot_worker_reduction_vs_slowdown(
    best_by_worker_count,
    reference_worker_count,
    maximum_slowdown_percent,
    output_folder,
):
    allowed_slowdowns = []
    worker_counts = []
    workers_removed_values = []
    actual_slowdowns = []

    reference_maximum = (
        best_by_worker_count[reference_worker_count]["maximum"]
    )

    for allowed_slowdown in range(
        1,
        maximum_slowdown_percent + 1,
    ):
        result = find_smallest_worker_count_near_reference(
            best_by_worker_count,
            reference_worker_count=reference_worker_count,
            allowed_slowdown_percent=float(
                allowed_slowdown
            ),
        )

        if result is None:
            continue

        workers_removed = (
            reference_worker_count
            - result["worker_count"]
        )

        actual_slowdown = (
            result["maximum"]
            / reference_maximum
            - 1.0
        ) * 100.0

        allowed_slowdowns.append(allowed_slowdown)
        worker_counts.append(result["worker_count"])
        workers_removed_values.append(workers_removed)
        actual_slowdowns.append(actual_slowdown)

    plt.figure(figsize=(14, 7))

    plt.plot(
        allowed_slowdowns,
        worker_counts,
        marker="o",
        linewidth=1.5,
    )

    plt.xlabel("Allowed slowdown relative to 127 workers (%)")
    plt.ylabel("Smallest acceptable worker count")
    plt.title(
        "Minimum Worker Count Versus Allowed Slowdown"
    )

    plt.xticks(range(1, maximum_slowdown_percent + 1))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "minimum_workers_vs_allowed_slowdown.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

   
# COST


def add_cost_metrics(
    benchmark_results,
    coordinator_count=1,
):
    """
    adds cost-related metrics to every worker-count/method result.

    definitions:

        total_cores = worker_count + coordinator_count

        multiplication_cost:
            total_cores * maximum_load

            This represents core-time cost if runtime is proportional
            to the largest worker load. Lower is better.

        division_score:
            total_cores / maximum_load

            This follows the proposed equation exactly. Higher is better,
            but it behaves more like an efficiency score than a cost.

        load_per_core:
            maximum_load / total_cores

            Lower is better.
    """

    for worker_count, worker_data in benchmark_results.items():
        total_cores = worker_count + coordinator_count

        for method_name, result in worker_data["methods"].items():
            maximum_load = result["maximum"]

            result["total_cores"] = total_cores

            result["multiplication_cost"] = (
                total_cores * maximum_load
            )

            if maximum_load > 0:
                result["division_score"] = (
                    total_cores / maximum_load
                )
            else:
                result["division_score"] = 0.0

            result["load_per_core"] = (
                maximum_load / total_cores
            )
        
def find_lowest_cost_method_by_worker_count(
    benchmark_results,
):
    """
    finds the method with the lowest core-time cost at each worker count.

    cost is defined as:

        (worker_count + 1) * largest_worker_load

    Lower is better.
    """

    cheapest_by_worker_count = {}

    for worker_count, worker_data in benchmark_results.items():
        best_method, best_result = min(
            worker_data["methods"].items(),
            key=lambda item: item[1]["multiplication_cost"],
            #key=lambda item: item[1]["division_score"],
        )

        cheapest_by_worker_count[worker_count] = {
            "method": best_method,
            "worker_count": worker_count,
            "total_cores": best_result["total_cores"],
            "maximum_load": best_result["maximum"],
            "cost": best_result["multiplication_cost"],
            #"cost": best_result["division_score"],
            "division_score": best_result["division_score"],
            "efficiency": best_result["efficiency"],
        }

    return cheapest_by_worker_count

def find_overall_lowest_cost_configuration(
    benchmark_results,
):
    """
    finds the lowest-cost combination of:

        worker count
        distribution method

    using:

        cost = (workers + 1) * maximum_load
    """

    best_configuration = None

    for worker_count, worker_data in benchmark_results.items():
        for method_name, result in worker_data["methods"].items():
            candidate = {
                "worker_count": worker_count,
                "total_cores": result["total_cores"],
                "method": method_name,
                "maximum_load": result["maximum"],
                "cost": result["multiplication_cost"],
                #"cost": result["division_score"],
                "division_score": result["division_score"],
            }

            if (
                best_configuration is None
                or candidate["cost"] < best_configuration["cost"]
            ):
                best_configuration = candidate

    return best_configuration


def save_cost_comparison_csv(
    comparisons,
    output_folder,
):
    output_path = os.path.join(
        output_folder,
        "runtime_and_cost_comparison.csv",
    )

    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "worker_count",
            "total_cores_including_coordinator",
            "workers_removed",
            "cores_removed",
            "best_method",
            "maximum_load",

            "estimated_runtime_percent",
            "estimated_runtime_change_percent",

            "core_load_cost",
            "cost_percent_of_reference",
            "cost_change_percent",
            "cost_savings_percent",

            "reference_method",
            "reference_maximum_load",
            "reference_total_cores",
            "reference_core_load_cost",
        ])

        for result in comparisons:
            writer.writerow([
                result["worker_count"],
                result["total_cores"],
                result["workers_removed"],
                result["cores_removed"],
                result["method"],
                result["maximum_load"],

                result["runtime_percent"],
                result["runtime_change_percent"],

                result["cost"],
                result["cost_percent"],
                result["cost_change_percent"],
                result["cost_savings_percent"],

                result["reference_method"],
                result["reference_maximum_load"],
                result["reference_total_cores"],
                result["reference_cost"],
            ])

    print(
        f"Runtime and cost comparison saved to: "
        f"{output_path}"
    )      
        


def print_cost_saving_options(
    comparisons,
    minimum_savings_percent=0.0,
    maximum_runtime_increase_percent=30.0,
):
    """
    prints configurations whose core-load cost is lower than the
    reference while remaining within the allowed runtime increase.
    """

    valid_options = [
        result
        for result in comparisons
        if (
            result["cost_savings_percent"]
            >= minimum_savings_percent
            and result["runtime_change_percent"]
            <= maximum_runtime_increase_percent
        )
    ]

    valid_options.sort(
        key=lambda result: result["worker_count"],
        reverse=True,
    )

    print()
    print("=" * 160)
    print("RUNTIME AND CORE-LOAD COST COMPARISON")
    print("=" * 160)

    heading = (
        f"{'Workers':>8}"
        f"{'Cores':>8}"
        f"{'Removed':>10}"
        f"{'Method':>20}"
        f"{'Runtime %':>14}"
        f"{'Runtime change':>18}"
        f"{'Cost %':>14}"
        f"{'Cost savings':>16}"
        f"{'Core-load cost':>20}"
        f"{'Maximum load':>18}"
    )

    print(heading)
    print("-" * len(heading))

    for result in valid_options:
        print(
            f"{result['worker_count']:>8}"
            f"{result['total_cores']:>8}"
            f"{result['cores_removed']:>10}"
            f"{result['method']:>20}"
            f"{result['runtime_percent']:>13.2f}%"
            f"{result['runtime_change_percent']:>17.2f}%"
            f"{result['cost_percent']:>13.2f}%"
            f"{result['cost_savings_percent']:>15.2f}%"
            f"{result['cost']:>20.2f}"
            f"{result['maximum_load']:>18.2f}"
        )

def compare_cost_to_reference(
    benchmark_results,
    reference_worker_count=127,
):
    """
    compares the best method at every worker count against the best
    method at the reference worker count.

    runtime:
        proportional to maximum worker load

    cost:
        (worker_count + 1 coordinator) * maximum worker load

    lower cost is better.
    """

    if reference_worker_count not in benchmark_results:
        raise ValueError(
            f"Reference worker count "
            f"{reference_worker_count} was not benchmarked."
        )

    reference_methods = (
        benchmark_results[reference_worker_count]["methods"]
    )

    # Select the fastest method at the reference worker count.
    reference_method, reference_result = min(
        reference_methods.items(),
        key=lambda item: item[1]["maximum"],
    )

    reference_maximum_load = reference_result["maximum"]
    reference_total_cores = reference_result["total_cores"]
    reference_cost = reference_result["multiplication_cost"]

    comparisons = []

    for worker_count, worker_data in benchmark_results.items():

        # Select the fastest method for this worker count.
        best_method, best_result = min(
            worker_data["methods"].items(),
            key=lambda item: item[1]["maximum"],
        )

        current_maximum_load = best_result["maximum"]
        current_total_cores = best_result["total_cores"]
        current_cost = best_result["multiplication_cost"]

        runtime_percent = (
            current_maximum_load
            / reference_maximum_load
            * 100.0
        )

        runtime_change_percent = (
            runtime_percent - 100.0
        )

        cost_percent = (
            current_cost
            / reference_cost
            * 100.0
        )

        cost_change_percent = (
            cost_percent - 100.0
        )

        cost_savings_percent = (
            100.0 - cost_percent
        )

        comparisons.append({
            "worker_count": worker_count,
            "total_cores": current_total_cores,

            "workers_removed": (
                reference_worker_count - worker_count
            ),

            "cores_removed": (
                reference_total_cores - current_total_cores
            ),

            "method": best_method,
            "maximum_load": current_maximum_load,

            "runtime_percent": runtime_percent,
            "runtime_change_percent": runtime_change_percent,

            "cost": current_cost,
            "cost_percent": cost_percent,
            "cost_change_percent": cost_change_percent,
            "cost_savings_percent": cost_savings_percent,

            "reference_method": reference_method,
            "reference_maximum_load": reference_maximum_load,
            "reference_total_cores": reference_total_cores,
            "reference_cost": reference_cost,
        })

    return comparisons

def find_best_configuration_for_savings(
    comparisons,
    required_cost_savings_percent,
    maximum_runtime_increase_percent,
):
    """
    finds the fastest configuration satisfying both:

        cost savings >= requested amount
        runtime increase <= permitted amount
    """

    valid_options = [
        result
        for result in comparisons
        if (
            result["cost_savings_percent"]
            >= required_cost_savings_percent
            and result["runtime_change_percent"]
            <= maximum_runtime_increase_percent
        )
    ]

    if not valid_options:
        return None

    return min(
        valid_options,
        key=lambda result: result["maximum_load"],
    )
    
def plot_cost_vs_worker_count(
    benchmark_results,
    output_folder,
):
    worker_counts = np.asarray(
        sorted(benchmark_results.keys()),
        dtype=int,
    )

    first_worker_data = benchmark_results[
        worker_counts[0]
    ]

    method_names = list(
        first_worker_data["methods"].keys()
    )

    plt.figure(figsize=(15, 8))

    for method_name in method_names:
        costs = np.asarray(
            [
                benchmark_results[worker_count]
                ["methods"][method_name]
                ["multiplication_cost"]
                #["division_score"]
                #division_score
                for worker_count in worker_counts
            ],
            dtype=float,
        )

        plt.plot(
            worker_counts,
            costs,
            linewidth=1.5,
            label=method_name,
        )

    plt.xlabel("Number of worker CPUs")
    plt.ylabel("(Workers + coordinator) × maximum load")
    plt.title(
        "Estimated Computational Cost Versus Worker Count\n"
        "Lower values indicate lower estimated core-time cost"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "cost_vs_worker_count.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

def plot_runtime_vs_cost_savings(
    comparisons,
    output_folder,
):
    runtime_changes = np.asarray(
        [
            result["runtime_change_percent"]
            for result in comparisons
        ],
        dtype=float,
    )

    cost_savings = np.asarray(
        [
            result["cost_savings_percent"]
            for result in comparisons
        ],
        dtype=float,
    )

    worker_counts = np.asarray(
        [
            result["worker_count"]
            for result in comparisons
        ],
        dtype=int,
    )

    plt.figure(figsize=(14, 8))

    plt.scatter(
        runtime_changes,
        cost_savings,
        s=25,
    )

    # Label selected points to avoid overcrowding.
    for index, worker_count in enumerate(worker_counts):
        if (
            worker_count % 5 == 0
            or worker_count == 127
            or worker_count == 1
        ):
            plt.annotate(
                str(worker_count),
                (
                    runtime_changes[index],
                    cost_savings[index],
                ),
                fontsize=8,
            )

    plt.axhline(
        0.0,
        linestyle="--",
        linewidth=1,
    )

    plt.axvline(
        0.0,
        linestyle="--",
        linewidth=1,
    )

    plt.xlabel("Predicted runtime increase versus 127 workers (%)")
    plt.ylabel("Estimated cost savings versus 127 workers (%)")

    plt.title(
        "Runtime–Cost Tradeoff\n"
        "Point labels show worker count"
    )

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = os.path.join(
        output_folder,
        "runtime_vs_cost_savings.png",
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# MAIN


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    
    (
        R_values,
        epsilon_values,
        weight_values,
        load_factors,
    ) = read_four_row_file(INPUT_FILE)

    print(f"Read input file: {INPUT_FILE}")
    print(f"Load factor at epsilon bin 0: {load_factors[0]}")
    #reads input data

    tasks = create_tasks(load_factors)
    #creates the collision integral task

    expected_task_count = NUMBER_OF_EPSILON_BINS * 2

    if len(tasks) != expected_task_count:
        raise RuntimeError(
            f"Expected {expected_task_count} collision-integral tasks, "
            f"but found {len(tasks)}."
        )

    print(f"Created {len(tasks)} collision-integral tasks.")

    selected_worker_count = NUMBER_OF_WORKERS

    print()
    print("=" * 100)
    print(
        f"DETAILED DISTRIBUTION COMPARISON: "
        f"{selected_worker_count} WORKERS"
    )
    print("=" * 100)

    (
        assignments_by_method,
        results,
        theoretical_lower_bound,
    ) = run_all_distribution_methods(
        tasks,
        selected_worker_count,
    )

    print_results(results, theoretical_lower_bound,)

    for method_name, assignments in assignments_by_method.items():
        safe_method_name = (
            method_name.lower()
            .replace("-", "_")
            .replace(" ", "_")
        )

        save_assignments_csv(
            safe_method_name,
            assignments,
            OUTPUT_FOLDER,
        )
    # save assignment CSV files for the selected worker count.

    save_summary_csv(
        results,
        OUTPUT_FOLDER,
        theoretical_lower_bound,
    )
    # save summary and graphs for the selected worker count.
    
    plot_worker_loads(
        results,
        OUTPUT_FOLDER,
    )

    plot_sorted_worker_loads(
        results,
        OUTPUT_FOLDER,
    )

    plot_maximum_loads(
        results,
        OUTPUT_FOLDER,
        theoretical_lower_bound,
    )

    plot_task_counts(
        results,
        assignments_by_method,
        OUTPUT_FOLDER,
    )

    minimum_worker_count = 1
    maximum_worker_count = NUMBER_OF_WORKERS
    # BENCHMARK ALL WORKER COUNTS FROM 1 THROUGH 127

    print()
    print("=" * 100)
    print(
        f"BENCHMARKING WORKER COUNTS "
        f"{minimum_worker_count} THROUGH "
        f"{maximum_worker_count}"
    )
    print("=" * 100)

    benchmark_results = benchmark_worker_counts(
        tasks,
        minimum_workers=minimum_worker_count,
        maximum_workers=maximum_worker_count,
    )
    '''
     ADD COST METRICS:

     adds:
    
     total_cores = workers + 1 coordinator
    
     multiplication_cost =
         total_cores * maximum worker load
    
     division_score =
         total_cores / maximum worker load
    
     the multiplication value is the more natural core-time
     cost estimate. The division value preserves the proposed
     equation for comparison.
    '''
    add_cost_metrics(
        benchmark_results,
        coordinator_count=1,
    )

    best_by_worker_count = (
        find_best_method_for_each_worker_count(
            benchmark_results
        )
    )
    # finds the best method for each worker count

    (
        method_win_counts,
        ties_by_worker_count,
    ) = count_method_wins(
        benchmark_results
    )
    
    # SAVE BENCHMARK RESULTS

    save_worker_count_benchmark_csv(
        benchmark_results,
        OUTPUT_FOLDER,
    )

    save_best_method_csv(
        best_by_worker_count,
        OUTPUT_FOLDER,
    )

    save_method_win_counts(
        method_win_counts,
        OUTPUT_FOLDER,
    )

    # CREATE BENCHMARK GRAPHS

    plot_maximum_load_vs_worker_count(
        benchmark_results,
        OUTPUT_FOLDER,
    )

    plot_relative_gap_vs_worker_count(
        benchmark_results,
        OUTPUT_FOLDER,
    )

    plot_efficiency_vs_worker_count(
        benchmark_results,
        OUTPUT_FOLDER,
    )

    plot_best_tested_load_vs_worker_count(
        best_by_worker_count,
        OUTPUT_FOLDER,
    )

    # PRINT METHOD WIN COUNTS

    print()
    print("=" * 100)
    print("METHOD WINS OR TIES ACROSS WORKER COUNTS")
    print("=" * 100)

    sorted_method_wins = sorted(
        method_win_counts.items(),
        key=lambda item: item[1],
        reverse=True,
    )

    for method_name, win_count in sorted_method_wins:
        print(
            f"{method_name:<22} "
            f"{win_count:>3} worker counts"
        )

    # REDUCED-WORKER COMPARISON FROM 1% THROUGH 30%

    reference_worker_count = NUMBER_OF_WORKERS
    minimum_slowdown_percent = 1
    maximum_slowdown_percent = 30

    print_reduced_worker_comparison(
        best_by_worker_count,
        reference_worker_count=reference_worker_count,
        minimum_slowdown_percent=minimum_slowdown_percent,
        maximum_slowdown_percent=maximum_slowdown_percent,
    )

    # SAVE REDUCED-WORKER RESULTS

    save_reduced_worker_comparison(
        best_by_worker_count,
        reference_worker_count=reference_worker_count,
        minimum_slowdown_percent=minimum_slowdown_percent,
        maximum_slowdown_percent=maximum_slowdown_percent,
        output_folder=OUTPUT_FOLDER,
    )

    plot_worker_reduction_vs_slowdown(
        best_by_worker_count,
        reference_worker_count=reference_worker_count,
        maximum_slowdown_percent=maximum_slowdown_percent,
        output_folder=OUTPUT_FOLDER,
    )

    # COST COMPARISON AGAINST 127 WORKERS

    cost_comparisons = compare_cost_to_reference(
        benchmark_results,
        reference_worker_count=reference_worker_count,
    )

    save_cost_comparison_csv(
        cost_comparisons,
        OUTPUT_FOLDER,
    )

    print_cost_saving_options(
        cost_comparisons,
        minimum_savings_percent=0.0,
        maximum_runtime_increase_percent=30.0,
    )

    # FIND THE OVERALL LOWEST-COST CONFIGURATION

    cheapest_by_worker_count = (
        find_lowest_cost_method_by_worker_count(
            benchmark_results
        )
    )

    overall_cheapest = (
        find_overall_lowest_cost_configuration(
            benchmark_results
        )
    )

    print()
    print("=" * 100)
    print("OVERALL LOWEST-COST CONFIGURATION")
    print("=" * 100)

    print(
        f"Workers:       "
        f"{overall_cheapest['worker_count']}"
    )

    print(
        f"Total cores:   "
        f"{overall_cheapest['total_cores']}"
    )

    print(
        f"Method:        "
        f"{overall_cheapest['method']}"
    )

    print(
        f"Maximum load:  "
        f"{overall_cheapest['maximum_load']:.2f}"
    )
    
    print(
        f"Core-load cost: "
        f"{overall_cheapest['cost']:.2f}"
    )

    print(
        f"Division score: "
        f"{overall_cheapest['division_score']:.12f}"
        #f"{overall_cheapest['multiplication_cost']:.12f}"
        #multiplication_cost
    )
    
    # EXAMPLE COST-SAVING TARGET

    required_cost_savings_percent = 10.0
    allowed_runtime_increase_percent = 30.0

    savings_option = find_best_configuration_for_savings(
        cost_comparisons,
        required_cost_savings_percent=(
            required_cost_savings_percent
        ),
        maximum_runtime_increase_percent=(
            allowed_runtime_increase_percent
        ),
    )

    print()
    print("=" * 100)
    print("TARGET COST-SAVING CONFIGURATION")
    print("=" * 100)

    if savings_option is None:
        print(
            f"No tested configuration saves at least "
            f"{required_cost_savings_percent:.2f}% "
            f"while staying within "
            f"{allowed_runtime_increase_percent:.2f}% "
            f"runtime increase."
        )
    else:
        print(
            f"Required cost savings: "
            f"{required_cost_savings_percent:.2f}%"
        )

        print(
            f"Maximum runtime increase: "
            f"{allowed_runtime_increase_percent:.2f}%"
        )

        print(
            f"Workers:       "
            f"{savings_option['worker_count']}"
        )

        print(
            f"Total cores:   "
            f"{savings_option['total_cores']}"
        )

        print(
            f"Method:        "
            f"{savings_option['method']}"
        )

        print(
            f"Runtime change: "
            f"{savings_option['runtime_change_percent']:.2f}%"
        )

        print(
            f"Cost savings:   "
            f"{savings_option['cost_savings_percent']:.2f}%"
        )

        print(
            f"Maximum load:   "
            f"{savings_option['maximum_load']:.2f}"
        )

    # CREATE COST GRAPHS

    plot_cost_vs_worker_count(
        benchmark_results,
        OUTPUT_FOLDER,
    )

    plot_runtime_vs_cost_savings(
        cost_comparisons,
        OUTPUT_FOLDER,
    )

    print()
    print("=" * 100)
    print("ANALYSIS COMPLETE")
    print("=" * 100)
    print(f"All results were saved in: {OUTPUT_FOLDER}")


if __name__ == "__main__":
    main()
    
'''
to give an understanding for each sorting method:

wall: 
    the wall method is nothing special, its the standard, sterotypical version of how people tend to 
    sort things. go from left to right, onces the limit is reached you start over, from the far left, 
    and work your way to the far right. much like how certain walls are built.

snake: 
    snake method isnt too different from the wall method. much like before you start from left to right. 
    the main difference occurs when you are at the far right. the moment you reach the far right, you change
    directions. not moving from left to right, to right to left. so itll look like this. left to right 1-127
    then at 127, that same core will also recieve workload 128. then core 126 will get 129, so on and so forth.
    the pattern then repeats like a snake slithering across the ground.
    
mountain:
    mountain method works far differently than the rest. essentially how this method works is by grabbing
    the theoretical largest summation values. say the last 2 workload. and gives back a number. now 
    cores 1 to 126 must relate to that workload. so if the sum of the last 2 workloads are 600k, then core 1
    will recieve a workload that is relatievely aroudn 600k, the moment it goes over a certain amount, it then
    gives the new workload to the next core, so on and so forth. 

largest first:
    largest first method isnt too complicated. basically the computer will sort out the workload from 
    largest to least then give the lowest workload core the next series of the workload array. so in 
    theory if we have a workload of the following, 1,3,2,5. itll assign the next work load to the core 
    with a workload of 1. then the next work load will be assigned to the core with the work load of 2.
    so on and so forth.

high-low snake:
    high low snake is a combination of the snake method and the largest first methods. essentially how this 
    sorting method will work, is by sorting the array from largest to smallest like before. however much like
    the snake method, it'll give the cores the workload from left to right, then right to left, and keep 
    alternating. this is to balance out the workload. if a core recieved the largest workload from the previous
    row, then this method will ensure it will recieve the smallest one in the next row. this is just to balance
    out the workload for each of the cores. 
    

the main 2 methods that are the most promising are largest first, and mountain. largest first method 
makes sense. literally just distributes the next largest workload to the next lowest workload on the core.
however, even that isnt perfect. it works most of the time really well. but the more cores or the more data 
(epsilons, weights, and r values) we have. the more and more the largest first method can potentially make.
meaning there will be gaps in which the mountain method will be better. meaning at points in time, the 
mountain method will work better than the largest first method. neither are reliable. but its possible to 
alternate between the 2, if things start to become a problem.    
    

'''
    