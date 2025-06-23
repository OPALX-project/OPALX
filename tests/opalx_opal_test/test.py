import subprocess
import time

import sys
sys.path.append('/psi/home/adelmann/git/pyOPALTools')
from opal import load_dataset

import numpy as np
import pandas as pd
import os

from template import *

import matplotlib.pyplot as plt

MULTI_THREAD_PREFIX = "multi_"

# Epsilon used by the tests
EPSILON = 1e-9

# beatiful colors
def color(r, g, b, background = False):
    return f"\033[38;2;{r};{g};{b}m"

RED = color(255, 0, 0)
MAGENTA = color(245, 66, 242)
BLUE = color(0, 0, 255)
GREEN = color(0, 255, 0)
ORANGE = color(255, 198, 77)
WHITE = "\033[0m"

BOLD="\033[1m"
RESET="\033[0m"

OPALX_EXECUTABLE_FILE = "/data/user/binder_j/opalx-elements/build_rel/src/opalx"
AMOUNT_THREADS = "16"

# define the parameters
parameters = {
    "run1" : {
        "amount" : "1e7",
        "avg" : 10,
        "ref" : "ref-33step.stat"
    }
}


# coloumns that get plotted in the end
plotting_cols = [
    "rms_x",
    "rms_y",
    "rms_s",
    "emit_x",
    "emit_y",
    "emit_s",
    "max_x",
    "max_y",
    "max_s",
    "xpx",
    "ypy",
    "zpz",
    "Dx",
    "DDx",
    "Dy",
    "DDy"
]


# init cols
for key in parameters.keys():
    parameters[key]["data"] = {}
    for cols in plotting_cols:
        parameters[key]["data"][cols] = {"error":[]}

def run_opalx(filename):
    if os.path.isfile(f"opalx/{filename.replace('in', 'stat')}"):
        print(f"{MAGENTA}Skipping execution with file {BLUE}{filename}{MAGENTA} for {ORANGE}Opalx{WHITE}")
        return

    try:
        print(f"Starting {ORANGE}Slurm script{WHITE} on dataset {BLUE + filename + WHITE}", end="", flush=True)

        result = subprocess.run(["sbatch", f"{filename.replace('in', 'slurm')}"], capture_output=True, text=True, check=True, cwd="opalx")

        if result.returncode == 0:
            print(f" --- {GREEN}Done{WHITE}")
        else:
            print(f"{RED}Return Code: {result.returncode + WHITE}")
            print(result.stdout)
            print(f"Code execution took {BOLD}{stop_time - start_time:.2f} seconds{RESET}")

    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}:")
        print(f"Standard Error:\n{e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Creates slurm and input file
def create_file(filename, amount, steps, seed):
    print(f"Creating {BLUE+filename+WHITE} and {BLUE+filename.replace('in', 'slurm')+WHITE}")
    print(f"Creating {BLUE}{MULTI_THREAD_PREFIX}{filename+WHITE} and {BLUE}{MULTI_THREAD_PREFIX}{filename.replace('in', 'slurm')+WHITE}")

    with open(f"opalx/{filename}", "w") as f:
        f.write(get_opalx_string(amount, steps, str(seed)))
        f.close()
    
    with open(f"opalx/{filename.replace('in', 'slurm')}", "w") as f:
        f.write(get_slurm_string(OPALX_EXECUTABLE_FILE, filename))
        f.close()

    with open(f"opalx/{MULTI_THREAD_PREFIX}{filename}", "w") as f:
        f.write(get_opalx_string(amount, steps, str(seed * 100)))
        f.close()

    with open(f"opalx/{MULTI_THREAD_PREFIX}{filename.replace('in', 'slurm')}", "w") as f:
        f.write(get_slurm_string(OPALX_EXECUTABLE_FILE, f"{MULTI_THREAD_PREFIX}{filename}", AMOUNT_THREADS))
        f.close()

def compare(stat1, stat2):
    amount = len(stat1.columns)
    passed = 0
    for c in stat2.columns:
        if all(abs(stat1[c] - stat2[c]) < EPSILON) == False and c != "numParticles":
            print(f"\033[31mFailed {c} test: Max diff {max(abs(stat1[c] - stat2[c])):e}\033[0m")
            pass
        else:
            #print(f"\033[34mPassed {c} test\033[0m")  
            passed += 1

    ratio = passed/amount
    
    print(f"\033[35m --- Passed {ratio:.2%} of tests --- \033[0m")


def merge_data(run_name):
    number_of_run = parameters[run_name]["avg"]

    # first load in the first data set as a reference for the cols
    filename = f"{key}-{param['amount']}-{0}.stat"
    data = load_dataset("opalx", fname=filename).dataframe
    # start from 1 instead of 0
    for i in range(1, number_of_run):
        filename = f"{key}-{param['amount']}-{i}.stat"
        data += load_dataset("opalx", fname=filename).dataframe
    print(f"Merged {BLUE}{key}-{param['amount']}-xxx.stat{WHITE}")

    # at the end divide by the amount of run
    data /= number_of_run

    # now do the same for the multithreaded one
    # first load in the first data set as a reference for the cols
    filename = f"{MULTI_THREAD_PREFIX}{key}-{param['amount']}-{0}.stat"
    multidata = load_dataset("opalx", fname=filename).dataframe
    # start from 1 instead of 0
    for i in range(1, number_of_run):
        filename = f"{MULTI_THREAD_PREFIX}{key}-{param['amount']}-{i}.stat"
        multidata += load_dataset("opalx", fname=filename).dataframe
    print(f"Merged {BLUE}{MULTI_THREAD_PREFIX}{key}-{param['amount']}-xxx.stat{WHITE}")

    # at the end divide by the amount of run
    multidata /= number_of_run

    return (data, multidata)


def compare_data(run_name):
    data, multi_data = merge_data(run_name)
    # load in reference
    reference_data = load_dataset("reference", fname=parameters[run_name]["ref"]).dataframe

    print(f"Comparing single threaded data")
    compare(data, reference_data)
    print(f"Comparing multi threaded data")
    compare(multi_data, reference_data)

    print(f"Plotting for run {BLUE}{run_name}{WHITE}")
    plot_all(data, multi_data, reference_data, run_name)

def plot_all(data, multi_data, reference, name):
    for col in plotting_cols:
        x = np.arange(1, len(data["t"]) + 1)

        plt.figure(figsize=(14, 7))
        plt.suptitle(f"Comparison to OPAL with {reference['numParticles'][0]:.0e} particles @ $\Delta t$ = {reference['t'][0]} ns\n OPALX simulation using {data['numParticles'][0]:.0e} particles averaged over {parameters[name]['avg']} runs")
        plt.subplot(121)

        plt.plot(x, data[col], label="Single threaded")
        plt.plot(x, multi_data[col], label=f"{AMOUNT_THREADS} cores")
        plt.plot(x, reference[col], ls="--", label="Reference")

        plt.grid()
        plt.legend()
        plt.xlabel("Steps")
        plt.title(f"Raw data from col {col}")

        plt.subplot(122)
        plt.semilogy(x, np.abs(data[col] - reference[col]), label="Single threded")
        plt.semilogy(x, np.abs(multi_data[col] - reference[col]), label=f"{AMOUNT_THREADS} cores")

        plt.title("Absolute error compared to the reference")
        plt.xlabel("Steps")
        plt.ylabel("Absolute Error")
        plt.grid()
        plt.legend()

        plt.savefig(f"{name}-{col}.png", dpi=300)
        plt.close()

def calculate_statistics():
    # get all runs
    for (key, param) in parameters.items():
        # for each run go over all observered data points
        for data_key in param["data"].keys():
            # now calculate any statisitcs like mean or median
            param["data"][data_key]["mean_error"]   = np.mean(param["data"][data_key]["error"], axis=1)
            param["data"][data_key]["median_error"] = np.median(param["data"][data_key]["error"], axis=1)

def watch(filename, steps):
    # get the filenames and repplace the *.in to *.stat
    opalx_filename = filename.replace("in", "stat")
    print(f"Waiting for {BLUE}{opalx_filename}{WHITE}", end="", flush=True)
    while True:
        if os.path.isfile(f"opalx/{opalx_filename}") == False:
            time.sleep(3)
            continue

        opalxstat = load_dataset('opalx', fname=opalx_filename).dataframe
        if len(opalxstat["t"]) == steps:
            break
        time.sleep(3)

    print(f" --- {GREEN}Found{WHITE}")


def get_user_choice():
    while True:
        user_input = input("Run using slurm, locally or cancel ([s],c): ").lower().strip()
        if user_input == 's' or user_input == '': # default
            print("Running using Slurm")
            return 'slurm'
        elif user_input == 'c':
            return "compare"
        else:
            print(f"{RED}Invalid choice.{RESET} Please enter 's' or 'c'.")

if __name__ == "__main__":
    choice = get_user_choice()

    if choice == "slurm":
        for (key, param) in parameters.items():
            # get the amount of steps from the reference
            opal_filename  = param["ref"]
            opalstat = load_dataset('reference', fname=opal_filename).dataframe
            steps = len(opalstat["t"])
            parameters[key]["steps"] = steps
            
            # 
            for run_number in range(param["avg"]):
                filename = f"{key}-{param['amount']}-{run_number}.in"
                create_file(filename, param['amount'], steps, run_number + 1)
                
                run_opalx(filename)
                run_opalx(f"{MULTI_THREAD_PREFIX}{filename}")

            for run_number in range(param["avg"]):
                filename = f"{key}-{param['amount']}-{run_number}.in"
                watch(filename, steps)
                watch(f"{MULTI_THREAD_PREFIX}{filename}", steps)

    for (key, param) in parameters.items():
        compare_data(key)
    