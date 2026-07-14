import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "./analysis")
sys.path.append("./analysis")
sys.path.insert(0, "../")
sys.path.append("../")
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

def debug_1_line(eps, eps2, data, data2, row_index, index2, title, filename, ylabel):
    plt.figure(figsize=(8, 6))

    plt.semilogy(eps, data[row_index, :], label="positive")
    plt.semilogy(eps, -data[row_index, :], label="negative")
    
    plt.semilogy(eps2, data2[index2, :], linestyle="--", label="original positive")
    plt.semilogy(eps2, -data2[index2, :], linestyle="--", label="original negative")

    plt.grid(True)
    plt.xlabel("epsilon")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc="best")

    plt.savefig(os.path.join(output_folder, filename), dpi=300, bbox_inches="tight")
    plt.close()
#this code above is used for the default cases, showing the original values unmodified before the step changes
#and plots them

def debug_6_lines(eps, eps2, data, data2, default_row, default_row2, row_indices, row_indices2, title, filename, ylabel):
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 12))
    axes = axes.flatten()

    for ax, row_index in zip(axes, row_indices):
        
        ax.semilogy(eps, data[row_index, :], label=f"modified row {row_index} positive")
        ax.semilogy(eps, -data[row_index, :], label=f"modified row {row_index} negative")
        # modified case
        
        ax.semilogy(eps, data[default_row, :], linestyle="--", label="default positive")
        ax.semilogy(eps, -data[default_row, :], linestyle="--", label="default negative")
        # default case overlay
        #   this results in a dashed line
        
    for ax, row_index in zip(axes, row_indices2):
        ax.semilogy(eps2, data2[row_index, :], linestyle="-.", label=f"modified og row {row_index} positive")
        ax.semilogy(eps2, -data2[row_index, :], linestyle="-.", label=f"modified og row {row_index} negative")
        # original modified case
        #   this results in a dash dot line
        
        ax.semilogy(eps2, data2[default_row2, :], linestyle=":", label="default og positive")
        ax.semilogy(eps2, -data2[default_row2, :], linestyle=":", label="default og negative")
        # default cases from the original data
        #   this results in a dotted line
        
        
        ax.grid(True)
        ax.set_xlabel("epsilon")
        ax.set_ylabel(ylabel)
        ax.set_title(f"row {row_index} vs default")
        ax.legend(loc="best", fontsize=7)

    fig.suptitle(title)
    fig.tight_layout()

    plt.savefig(os.path.join(output_folder, filename), dpi=300, bbox_inches="tight")
    plt.close()
#this code above, is mostly just used for plotting 
#   all different 6 cases in the test file

def debug_all_line(eps, eps2, data, data2, data_name):
    #rows_per_dataset = 5
    rows_per_dataset = 5
    number_of_datasets = 5
    
    for dataset_index in range(number_of_datasets):
        start = 4 + dataset_index * rows_per_dataset
        start2 = 2 + dataset_index * rows_per_dataset
        
        default_row = start
        default_row2 = start2
        
        modified_group_1 = range(start + 1, start + 6)
        modified_group_2 = range(start2 + 1, start2 + 6)

        dataset_number = dataset_index + 1

        debug_1_line(
            eps,
            eps2,
            data,
            data2,
            default_row,
            default_row2,
            f"Dataset {dataset_number} Default {data_name} vs Epsilon",
            f"dataset_{dataset_number}_default_{data_name}.png",
            data_name
        )

        debug_6_lines(
            eps,
            eps2,
            data,
            data2,
            default_row,
            default_row2,
            modified_group_1,
            modified_group_2,
            f"Dataset {dataset_number} Modified Group 1 {data_name} vs Default",
            f"dataset_{dataset_number}_modified_group_1_{data_name}.png",
            data_name
        )
#this code above plots all of the graphs so it can be shown

def debug_reform(data_file):
    results = {}

    eps_data = np.loadtxt(
        data_file,
        delimiter=",",
        max_rows=2,
        usecols=range(0, 206)
        #change that 2nd number in range, depending on the 
        #   number of epsilon values             
    )

    eps = eps_data[0]
    weights = eps_data[1]

    num_bins = len(eps)

    data = np.genfromtxt(
        data_file,
        delimiter=",",
        skip_header=2
    )
        #grabs the data from the data file

    data = data[~np.isnan(data).all(axis=1)]

    data = data[:, 2:]

    results["eps"] = eps
    results["w"] = weights
    results["N_bins"] = num_bins

    rho = data[:, 0:4 * num_bins]

    results["rho"] = rho

    return results
#this code is here to help overlay the original csv file over
#   the modified version, however the original had 206 epsilon values
#   the modified has about a 1k epsilone values

def plot_six_lines(eps, data, default_row, row_indices, title, filename, ylabel):
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 12))
    axes = axes.flatten()

    for ax, row_index in zip(axes, row_indices):
        
        ax.semilogy(eps, data[row_index, :], label=f"modified row {row_index} positive")
        ax.semilogy(eps, -data[row_index, :], label=f"modified row {row_index} negative")
        # modified case
        
        ax.semilogy(eps, data[default_row, :], linestyle="--", label="default positive")
        ax.semilogy(eps, -data[default_row, :], linestyle="--", label="default negative")
        # default case overlay

        ax.grid(True)
        ax.set_xlabel("epsilon")
        ax.set_ylabel(ylabel)
        ax.set_title(f"row {row_index} vs default")
        ax.legend(loc="best", fontsize=7)

    fig.suptitle(title)
    fig.tight_layout()

    plt.savefig(os.path.join(output_folder, filename), dpi=300, bbox_inches="tight")
    plt.close()
#this code above, is mostly just used for plotting 
#   all different 6 cases in the test file
       
def plot_all_datasets(eps, data, data_name):
    rows_per_dataset = 5
    number_of_datasets = 5

    for dataset_index in range(number_of_datasets):
        start = 4 + dataset_index * 5

        default_row = start
        modified_group_1 = range(start + 1, start + 6)

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
#this code above plots all of the graphs so it can be shown


def reform(data_file):
    results = {}

    eps_data = np.loadtxt(
        data_file,
        delimiter=",",
        max_rows=2,
        usecols=range(0, 1000)
        #change that 2nd number in range, depending on the 
        #   number of epsilon values             
    )

    eps = eps_data[0]
    weights = eps_data[1]

    num_bins = len(eps)

    data = np.genfromtxt(
        data_file,
        delimiter=",",
        skip_header=2
    )
        #grabs the data from the data file

    data = data[~np.isnan(data).all(axis=1)]

    data = data[:, 2:]

    results["eps"] = eps
    results["w"] = weights
    results["N_bins"] = num_bins

    rho = data[:, 0:4 * num_bins]

    results["rho"] = rho

    return results

#reform is a changed variation from the code 
#   process_data.py, because it does not know how 
#   graph something on the same csv file that has 
#   the epsilon, weights, and the density values all in 1 csv file

        
def make_P(rho):
    return rho[:,::4], rho[:,1::4], rho[:,2::4], rho[:,3::4]
#this code was taken from process_data.py 
#   was thinking about changing it, 
#   however it didnt seem like i needed to

def processes(extrapolation, original):
    os.makedirs(output_folder, exist_ok=True)

    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)

        if os.path.isfile(file_path):
            os.remove(file_path)
            
    data_file = "../" + extrapolation
    data_file_2= "../" + original
    data = reform(data_file)
    data2 = debug_reform(data_file_2)

    P0, Px, Py, Pz = make_P(data['rho'])
    P02, Px2, Py2, Pz2 = make_P(data2['rho'])

    eps_data = np.loadtxt(
                    data_file,
                    delimiter=',',
                    skiprows=0,
                    usecols=range(0, 1000) )
                    #change that 2nd number in range, depending on the 
                    #   number of epsilon values 
                    #   this one is using data_file as the epsilon values
                    
    eps_data_2 = np.loadtxt(
                    data_file_2,
                    delimiter=',',
                    skiprows=0,
                    usecols=range(0, 206) )
                    #this one is using the epsilon values from
                    #   data_file_2 so we can ensure the correct values
        
    eps = eps_data[0]
    eps2 = eps_data_2[0]
        
    debug_all_line(eps, eps2, P0, P02, "P0")
    debug_all_line(eps, eps2, Px, Px2, "Px")
    debug_all_line(eps, eps2, Py, Py2, "Py")
    debug_all_line(eps, eps2, Pz, Pz2, "Pz")
##############            
#this code above, is here just so that other files can run the code

RKextrap = sys.argv[1]
RKoriginal = sys.argv[2]
output_folder = "graphs_find_errors"
data_file = "../" + RKextrap
#data file that contains all the results of a extrapiolation test
data_file_2= "../" + RKoriginal
#data file that contains all the results of the original test
output_folder = "./solvingProblems/graphs_find_errors"
#just the folder where the graphs go

if os.path.exists(output_folder):
    shutil.rmtree(output_folder)

os.makedirs(output_folder)

    
os.makedirs(output_folder, exist_ok=True)

for filename in os.listdir(output_folder):
    file_path = os.path.join(output_folder, filename)

    if os.path.isfile(file_path):
        os.remove(file_path)    

#this is here to delete anything in the designated location
#just so we dont have an overflow of files
#remember to copy and paste the folder if we need any of the graphs

data = reform(data_file)
data2 = debug_reform(data_file_2)

P0, Px, Py, Pz = make_P(data['rho'])
P02, Px2, Py2, Pz2 = make_P(data2['rho'])

eps_data = np.loadtxt(
                data_file,
                delimiter=',',
                skiprows=0,
                usecols=range(0, 1000) )
                #change that 2nd number in range, depending on the 
                #   number of epsilon values 
                #   this one is using data_file as the epsilon values
                
eps_data_2 = np.loadtxt(
                data_file_2,
                delimiter=',',
                skiprows=0,
                usecols=range(0, 206) )
                #this one is using the epsilon values from
                #   data_file_2 so we can ensure the correct values


ar = [0, 1, 2, 3, 4, 5, 6]
binnumber = 1
#this may be unessecary code left over from a previous attempt
#but may need for debugging later

eps = eps_data[0]
eps2 = eps_data_2[0]
'''
#print(data)
print(eps)
print(len(eps))
print(len(P0[0]))
print(len(Px[0]))
print(len(Py[0]))
print(len(Pz[0]))
#   debugging statements to ensure 
#   everything is correct

print("original graphs\n")
print(eps2)
print(len(eps2))
print(len(P02[0]))
print(len(Px2[0]))
print(len(Py2[0]))
print(len(Pz2[0]))

#print(P02)
'''
#keep these debugging statements in case you need it
 
print("printing graphs, this may take a bit ...")
debug_all_line(eps, eps2, P0, P02, "P0")
#debug_all_line(eps, eps2, Px, Px2, "Px")
#debug_all_line(eps, eps2, Py, Py2, "Py")
#debug_all_line(eps, eps2, Pz, Pz2, "Pz")

#debugging graphs to overlay the original
#   and to overlay the new modified graphs

print("Graphs saved in:", output_folder)
