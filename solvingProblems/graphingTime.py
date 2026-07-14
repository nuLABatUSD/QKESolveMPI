import os
import re
import numpy as np
import matplotlib.pyplot as plt
import sys
#sys.path.insert(0, "C:\\Users\\astar\\Desktop\\workspace\\QKESolveMPI-main\\analysis")
sys.path.insert(0, "../analysis")
sys.path.append("../analysis")
#sys.path.insert(0, "C:/Users/astar/Desktop/workspace/QKESolveMPI-main/analysis")
#sys.path.append("C:/Users/astar/Desktop/workspace/QKESolveMPI-main/analysis")
#import process_data

import os
import shutil


def plot_one_line(eps, data, row_index, title, filename, ylabel):
    plt.figure(figsize=(8, 6))

    plt.semilogy(eps, data[row_index, :], label="positive")
    plt.semilogy(eps, -data[row_index, :], label="negative")

    plt.grid(True)
    plt.xlabel("epsilon")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc="best")

    plt.savefig(os.path.join(output_folder, filename), dpi=300, bbox_inches="tight")
    plt.close()
#this code above is used for the default cases, showing the original values unmodified before the step changes
#and plots them


def plot_six_lines(eps, data, default_row, row_indices, title, filename, ylabel):
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12))
    axes = axes.flatten()

    for ax, row_index in zip(axes, row_indices):
        # modified case
        ax.semilogy(eps, data[row_index, :], label=f"modified row {row_index} positive")
        ax.semilogy(eps, -data[row_index, :], label=f"modified row {row_index} negative")

        # default case overlay
        ax.semilogy(eps, data[default_row, :], linestyle="--", label="default positive")
        ax.semilogy(eps, -data[default_row, :], linestyle="--", label="default negative")

        ax.grid(True)
        ax.set_xlabel("epsilon")
        ax.set_ylabel(ylabel)
        ax.set_title(f"row {row_index} vs default")
        ax.legend(loc="best", fontsize=7)

    fig.suptitle(title)
    fig.tight_layout()

    plt.savefig(os.path.join(output_folder, filename), dpi=300, bbox_inches="tight")
    plt.close()
#this code above, is mostly just used for plotting all different 6 cases in the test file

def plot_all_datasets(eps, data, data_name):
    rows_per_dataset = 13
    number_of_datasets = 5

    for dataset_index in range(number_of_datasets):
        start = dataset_index * rows_per_dataset

        default_row = start
        modified_group_1 = range(start + 1, start + 7)
        modified_group_2 = range(start + 7, start + 13)

        dataset_number = dataset_index + 1

        plot_one_line(
            eps,
            data,
            default_row,
            f"Dataset {dataset_number} Default {data_name} vs Epsilon",
            f"dataset_{dataset_number}_default_{data_name}.png",
            data_name
        )

        plot_six_lines(
            eps,
            data,
            default_row,
            modified_group_1,
            f"Dataset {dataset_number} Modified Group 1 {data_name} vs Default",
            f"dataset_{dataset_number}_modified_group_1_{data_name}.png",
            data_name
        )

        plot_six_lines(
            eps,
            data,
            default_row,
            modified_group_2,
            f"Dataset {dataset_number} Modified Group 2 {data_name} vs Default",
            f"dataset_{dataset_number}_modified_group_2_{data_name}.png",
            data_name
        )
#this code above plots all of the graphs so it can be shown
def make_data_dictionary(data_file, eps_file):
    results = dict()
    
    data = np.loadtxt(data_file, delimiter=',')
    eps_data = np.loadtxt(eps_file, delimiter=',')
    num_bins = len(eps_data[0])

    results['N_bins'] = num_bins

    rho = data[:,2:2+4*num_bins]
    rhobar = data[:,2+num_bins*4:-2]
    results['rho'] = rho
    results['rhobar'] = rhobar

    P0 = rho[:,::4]
    Pz = rho[:,3::4]
    P0bar = rhobar[:,::4]
    Pzbar = rhobar[:,3::4]

    f_e = 0.5*P0*(1+Pz)
    f_m = 0.5*P0*(1-Pz)
    f_ebar = 0.5*P0bar*(1+Pzbar)
    f_mbar = 0.5*P0bar*(1-Pzbar)
    

    results['eps'] = eps_data[0]
    results['w'] = eps_data[1]

    dnde = []

    return results

def make_P(rho):
    return rho[:,::4], rho[:,1::4], rho[:,2::4], rho[:,3::4]


RK = sys.argv[1]
epsvalues = sys.argv[2]


output_folder = "graphs_output"

data_file = "./solvingProblems/" + RK
#data file that contains all the results of a test

epsilon_raw_values = "./solvingProblems/" + epsvalues
#the raw epsilon values from a data file

output_folder = "./solvingProblems/graphs_output"
#just the folder where the graphs go


if os.path.exists(output_folder):
    shutil.rmtree(output_folder)

os.makedirs(output_folder)

os.makedirs(output_folder, exist_ok=True)

for filename in os.listdir(output_folder):
    file_path = os.path.join(output_folder, filename)

    if os.path.isfile(file_path):
        os.remove(file_path)
#delete everything in the folder. so we know the code is running properly


data = make_data_dictionary(data_file, epsilon_raw_values)
P0, Px, Py, Pz = make_P(data['rho'])
eps_data = np.loadtxt(epsilon_raw_values, delimiter=',')

ar = [0, 1, 2, 3, 4, 5, 6]
binnumber = 1

eps = eps_data[0]

print("creating graphs, this may take a bit...")
#plot_all_datasets(eps, P0, "P0")
#plot_all_datasets(eps, Px, "Px")
#plot_all_datasets(eps, Py, "Py")
#plot_all_datasets(eps, Pz, "Pz")

#this code here is where all the inputs happen to create the graphs with the initial 
#   epsilon values 

print("Graphs saved in:", output_folder)


