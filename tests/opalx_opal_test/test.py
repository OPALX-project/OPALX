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

OPALX_EXECUTABLE_FILE = "/data/user/binder_j/opalx-elements/build/src/opalx"

# define the parameters
parameters = {
    "run1" : {
        "amount" : ["1e2", "1e3", "1e4", "1e5"],
        "ref" : "ref-10steps.stat"
    },
    "run2" : {
        "amount" : ["1e2", "1e3", "1e4", "5e4"],
        "ref" : "ref-33step.stat"
    }
}


# coloumns that get plotted in the end
plotting_cols = [
    "rms_x",
    "rms_y",
    "rms_s",
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
def create_file(filename, amount, steps):
    print(f"Creating {BLUE+filename+WHITE} and {BLUE+filename.replace('in', 'slurm')+WHITE}")

    with open(f"opalx/{filename}", "w") as f:
        f.write(get_opalx_string(amount, steps))
        f.close()
    
    with open(f"opalx/{filename.replace('in', 'slurm')}", "w") as f:
        f.write(get_slurm_string(OPALX_EXECUTABLE_FILE, filename))
        f.close()

def compare(stat1, stat2, runname):
    amount = len(stat1.columns)
    passed = 0
    for c in stat2.columns:
        if c in parameters[runname]["data"].keys():
            # add the error if its observedd
            parameters[runname]["data"][c]["error"].append((abs(stat1[c] - stat2[c])).to_numpy())

        if all(abs(stat1[c] - stat2[c]) < EPSILON) == False and c != "numParticles":
            print(f"\033[31mFailed {c} test: Max diff {max(abs(stat1[c] - stat2[c])):e}\033[0m")
            pass
        else:
            #print(f"\033[34mPassed {c} test\033[0m")  
            passed += 1

    ratio = passed/amount
    
    print(f"\033[35m --- Passed {ratio:.2%} of tests --- \033[0m")

def compare_data(filename, runname):
    # get the filenames and repplace the *.in to *.stat
    opalx_filename = filename.replace("in", "stat")
    opal_filename  = parameters[runname]["ref"]
    print(f"Loading stats {BLUE}opalx/{opalx_filename + WHITE} and {BLUE}reference/{opal_filename + WHITE}")
    # load in the dataset
    opalxstat = load_dataset('opalx', fname=opalx_filename).dataframe
    opalstat = load_dataset('reference', fname=opal_filename).dataframe
    print(f"Reference is using {MAGENTA}{opalstat['numParticles'][0]:.1e}{WHITE} Particles")
    # compare and save the errors
    compare(opalxstat, opalstat, runname)

def calculate_statistics():
    # get all runs
    for (key, param) in parameters.items():
        # for each run go over all observered data points
        for data_key in param["data"].keys():
            # now calculate any statisitcs like mean or median
            param["data"][data_key]["mean_error"]   = np.mean(param["data"][data_key]["error"], axis=1)
            param["data"][data_key]["median_error"] = np.median(param["data"][data_key]["error"], axis=1)

def plot_all():
    # get all runs
    for (run_name, param) in parameters.items():
        # load in the opal stats
        opal_filename  = param["ref"]
        opalstat = load_dataset('reference', fname=opal_filename).dataframe
        
        # for each run go over all obersaverbals
        for data_key in param["data"].keys():
            plt.figure(figsize=(15,5))
            plt.suptitle(f"Run '{run_name}'. Ref solution using {opalstat['numParticles'][0]:.0e} particles. $\Delta t$ = {opalstat['t'][0]} ns")
            x = np.arange(1, len(opalstat["t"]) + 1)
            plt.subplot(131)
            
            # iterate over all amounts again
            for amount in param["amount"]:
                # get the file name and stat
                opalx_filename = f"{run_name}-{amount}.stat"
                data = load_dataset('opalx', fname=opalx_filename).dataframe
                plt.plot(x, data[data_key], label=f"{amount}")
            plt.plot(x, opalstat[data_key], ls="--", label="Reference solution")
            plt.title(f"Raw {data_key} data")
            plt.xlabel("Step")
            plt.grid()
            plt.legend()


            plt.subplot(132)
            for (i, err) in enumerate(param["data"][data_key]["error"]):
                plt.semilogy(np.arange(1, len(err) + 1), err, label=f"{param['amount'][i]}")
            plt.title("Error semi log y")
            plt.xlabel("Steps")
            plt.ylabel("Error")
            plt.legend()
            plt.grid()

            plt.subplot(133)
            plt.title(f"Mean and median {data_key}")
            plt.loglog([float(s) for s in param["amount"]], param["data"][data_key]["mean_error"], label="Mean error")
            plt.loglog([float(s) for s in param["amount"]], param["data"][data_key]["median_error"], label="Median error")
            plt.xlabel("Particle amount")
            plt.ylabel("Error")
            plt.legend()
            plt.grid()

            plt.savefig(f"{run_name}-{data_key}.png")
            plt.close()

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
            
            # iterate over all particle sizes and create input and slurm file
            for amounts in param["amount"]:
                filename = f"{key}-{amounts}.in"
                create_file(filename, amounts, steps)
                run_opalx(filename)
            
            for amounts in param["amount"]:
                filename = f"{key}-{amounts}.in"
                watch(filename, steps)

    for (key, param) in parameters.items():
        for amounts in param["amount"]:
            filename = f"{key}-{amounts}.in"
            compare_data(filename, key)
    
    print("Plotting data...")
    calculate_statistics()
    plot_all()