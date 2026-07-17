import csv
import os
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# SETTINGS
# ============================================================

filename = "R_values-201-20-5.csv"
input_file = "./solvingProblems/" + filename

output_folder = "./solvingProblems/workerFolder"
output_filename = "current_worker_loads.png"

number_of_workers = 127
number_of_epsilon_bins = 205

# Two collision integrals per nonzero epsilon bin:
# one neutrino and one antineutrino.
number_of_collision_integrals = number_of_epsilon_bins * 2


def read_data_file(filename):
    """
    Reads a file containing four rows:

        Row 0: R
        Row 1: epsilon
        Row 2: quadrature weights
        Row 3: load factors

    The file may be comma-separated or whitespace-separated.
    """

    rows = []

    with open(filename, "r", newline="") as file:
        for raw_line in file:
            raw_line = raw_line.strip()

            if not raw_line:
                continue

            # Support either CSV or whitespace-separated data.
            if "," in raw_line:
                values = next(csv.reader([raw_line]))
            else:
                values = raw_line.split()

            cleaned_values = [
                float(value.strip())
                for value in values
                if value.strip() != ""
            ]

            rows.append(cleaned_values)

    if len(rows) != 4:
        raise ValueError(
            f"Expected exactly 4 nonempty rows, but found {len(rows)}."
        )

    row_lengths = [len(row) for row in rows]

    if any(length != 206 for length in row_lengths):
        raise ValueError(
            "Expected 206 values in every row.\n"
            f"Found row lengths: {row_lengths}"
        )

    R_values = np.asarray(rows[0], dtype=float)
    epsilon_values = np.asarray(rows[1], dtype=float)
    weight_values = np.asarray(rows[2], dtype=float)
    load_factors = np.asarray(rows[3], dtype=float)

    return R_values, epsilon_values, weight_values, load_factors


# ============================================================
# CURRENT WORK DISTRIBUTION
# ============================================================

def calculate_current_worker_loads(
    load_factors,
    worker_count,
    epsilon_bin_count
):
    """
    Reproduces the current supercomputer assignment.

    Global collision-integral tasks are numbered 1 through 410.

    Tasks 1-205:
        neutrino epsilon bins 1-205

    Tasks 206-410:
        antineutrino epsilon bins 1-205

    Worker w receives:

        w, w + worker_count, w + 2*worker_count, ...

    until the task number exceeds 410.
    """

    total_tasks = epsilon_bin_count * 2

    worker_loads = np.zeros(worker_count, dtype=float)
    worker_bins = []
    worker_tasks = []

    for worker in range(1, worker_count + 1):
        assigned_bins = []
        assigned_tasks = []
        total_load = 0.0

        task = worker

        while task <= total_tasks:
            # Convert global task number 1-410 into epsilon bin 1-205.
            epsilon_bin = ((task - 1) % epsilon_bin_count) + 1

            # load_factors[0] is epsilon bin 0.
            bin_load = load_factors[epsilon_bin]

            total_load += bin_load
            assigned_bins.append(epsilon_bin)
            assigned_tasks.append(task)

            task += worker_count

        worker_loads[worker - 1] = total_load
        worker_bins.append(assigned_bins)
        worker_tasks.append(assigned_tasks)

    return worker_loads, worker_bins, worker_tasks


# ============================================================
# SAVE RESULTS
# ============================================================

def save_worker_table(
    worker_loads,
    worker_bins,
    worker_tasks,
    output_path
):
    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "worker",
            "number_of_integrals",
            "task_numbers",
            "epsilon_bins",
            "total_load_factor"
        ])

        for index in range(len(worker_loads)):
            writer.writerow([
                index + 1,
                len(worker_tasks[index]),
                " ".join(map(str, worker_tasks[index])),
                " ".join(map(str, worker_bins[index])),
                worker_loads[index]
            ])


def plot_worker_loads(worker_loads, output_path):
    workers = np.arange(1, len(worker_loads) + 1)

    largest_index = int(np.argmax(worker_loads))
    smallest_index = int(np.argmin(worker_loads))

    largest_worker = largest_index + 1
    smallest_worker = smallest_index + 1

    largest_load = worker_loads[largest_index]
    smallest_load = worker_loads[smallest_index]
    average_load = np.mean(worker_loads)

    plt.figure(figsize=(14, 7))

    plt.plot(
        workers,
        worker_loads,
        marker="o",
        markersize=3,
        linewidth=1
    )

    plt.axhline(
        average_load,
        linestyle="--",
        linewidth=1,
        label=f"Average load = {average_load:.2f}"
    )

    plt.scatter(
        largest_worker,
        largest_load,
        s=70,
        zorder=3,
        label=(
            f"Largest: worker {largest_worker}, "
            f"load {largest_load:.2f}"
        )
    )

    plt.scatter(
        smallest_worker,
        smallest_load,
        s=70,
        zorder=3,
        label=(
            f"Smallest: worker {smallest_worker}, "
            f"load {smallest_load:.2f}"
        )
    )

    plt.xlabel("Worker number")
    plt.ylabel("Total load factor")
    plt.title(
        "Current Collision-Integral Load Distribution\n"
        f"{len(worker_loads)} workers, "
        f"{number_of_collision_integrals} collision integrals"
    )

    plt.xticks(np.arange(1, len(worker_loads) + 1, 5))
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    os.makedirs(output_folder, exist_ok=True)

    (
        R_values,
        epsilon_values,
        weight_values,
        load_factors
    ) = read_data_file(input_file)

    # The epsilon = 0 collision integral should not be calculated.
    print(f"Epsilon-bin-0 load factor: {load_factors[0]}")

    (
        worker_loads,
        worker_bins,
        worker_tasks
    ) = calculate_current_worker_loads(
        load_factors,
        number_of_workers,
        number_of_epsilon_bins
    )

    graph_path = os.path.join(output_folder, output_filename)
    table_path = os.path.join(
        output_folder,
        "current_worker_loads.csv"
    )

    plot_worker_loads(worker_loads, graph_path)

    save_worker_table(
        worker_loads,
        worker_bins,
        worker_tasks,
        table_path
    )

    largest_index = int(np.argmax(worker_loads))
    smallest_index = int(np.argmin(worker_loads))

    print("\nCurrent workload distribution")
    print("----------------------------------------")

    print(
        f"Largest load: {worker_loads[largest_index]}\n"
        f"Worker: {largest_index + 1}\n"
        f"Tasks: {worker_tasks[largest_index]}\n"
        f"Epsilon bins: {worker_bins[largest_index]}"
    )

    print()

    print(
        f"Smallest load: {worker_loads[smallest_index]}\n"
        f"Worker: {smallest_index + 1}\n"
        f"Tasks: {worker_tasks[smallest_index]}\n"
        f"Epsilon bins: {worker_bins[smallest_index]}"
    )

    print()
    print(f"Average load: {np.mean(worker_loads)}")
    print(f"Load range: {np.ptp(worker_loads)}")
    print(f"Standard deviation: {np.std(worker_loads)}")

    print()
    print(f"Graph saved to: {graph_path}")
    print(f"Worker table saved to: {table_path}")

    # Verification of the example in the problem statement.
    print()
    print("Worker 1 verification:")
    print(f"Tasks: {worker_tasks[0]}")
    print(f"Epsilon bins: {worker_bins[0]}")
    print(f"Total load: {worker_loads[0]}")


if __name__ == "__main__":
    main()


'''
4 rows
row 1 = r values
row 2 = epsilon
row 3 = weights 
row 4 = load factors

which is best to have a graph of first? perhaps load factor on the y axis and epsilon on x? 
'''
