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

OPALX_EXECUTABLE_FILE = "/data/user/binder_j/opalx-test/build/src/opalx"
MODE = "compare"

# define the parameters
parameters = {
    "run1" : {
        "steps" : "10",
        "amount" : ["1e2", "1e3", "1e4"],
        "ref" : "ref-10steps.stat"
    },
    #"run2" : {
    #    "steps" : "30",
    #    "amount" : ["1e2", "1e3", "1e4"],
    #    "ref" : "ref.stat"
    #}
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

    if MODE == "slurm":
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
    else:
        try:
            print(f"Running {ORANGE}Opalx{WHITE} on dataset {BLUE + filename + WHITE}", end="", flush=True)

            start_time = time.time()
            result = subprocess.run([OPALX_EXECUTABLE_FILE, f"{filename}", "--info", "10"], capture_output=True, text=True, check=True, cwd="opalx")
            stop_time = time.time()

            if result.returncode == 0:
                print(f" --- {GREEN}Done{WHITE} in {BOLD}{stop_time - start_time:.2f} seconds{RESET}")
            else:
                print(f"{RED}Return Code: {result.returncode + WHITE}")
                print(result.stdout)
                print(f"Code execution took {BOLD}{stop_time - start_time:.2f} seconds{RESET}")

        except subprocess.CalledProcessError as e:
            print(f"Command failed with return code {e.returncode}:")
            print(f"Standard Error:\n{e.stderr}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
        
def create_file(filename, amount, steps):
    if MODE == "slurm":
        print(f"Creating {BLUE+filename+WHITE} and {BLUE+filename.replace('in', 'slurm')+WHITE}")
    else:
        print(f"Creating {BLUE+filename+WHITE}")

    with open(f"opalx/{filename}", "w") as f:
        f.write(get_opalx_string(amount, steps))
        f.close()
    
    if MODE == "slurm":
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
    for (key, param) in parameters.items():
        # for each run go over all obersaverbals
        for data_key in param["data"].keys():
            plt.figure(figsize=(10,5))
            
            plt.subplot(121)
            for (i, err) in enumerate(param["data"][data_key]["error"]):
                plt.semilogy(np.arange(1, len(err) + 1), err, label=f"{param['amount'][i]}")
            plt.xlabel("Steps")
            plt.ylabel("Error")
            plt.title(f"{key} showing {data_key}")
            plt.legend()
            plt.grid()

            plt.subplot(122)
            plt.title(f"Mean and median {data_key}")
            plt.loglog([float(s) for s in param["amount"]], param["data"][data_key]["mean_error"], label="Mean error")
            plt.loglog([float(s) for s in param["amount"]], param["data"][data_key]["median_error"], label="Median error")
            plt.xlabel("Particle amount")
            plt.ylabel("Error")
            plt.legend()
            plt.grid()

            plt.savefig(f"{key}-{data_key}.png")
            plt.close()

def get_user_choice():
    global MODE
    while True:
        user_input = input("Run using slurm, locally or cancel (s,[l],c): ").lower().strip()
        if user_input == 's':
            MODE = "slurm"
            print("Running using Slurm")
            return 'slurm'
        elif user_input == 'l' or user_input == '': # 'l' or empty input (default)
            MODE = "local"
            print("Running on login node")
            return 'local'
        elif user_input == 'c':
            MODE = "compare"
            return "compare"
        else:
            print("Invalid choice. Please enter 's', 'l', or 'c'.")

def check_if_file_exists(filepath, pollingrate = 3, timeout = 60):
    start_time = time.time()
    print(f"Waiting for file '{BLUE+filepath+WHITE}'", end="", flush=True)
    while time.time() - start_time < timeout:
        if os.path.exists(filepath):
            print(f" --- {GREEN} Found {WHITE}")
            return True
        time.sleep(pollingrate)
    print(f" --- {RED} Timeout waiting for file '{filepath}' {WHITE}.")
    return False


if __name__ == "__main__":
    get_user_choice()

    if MODE == "slurm":
        for (key, param) in parameters.items():
            for amounts in param["amount"]:
                filename = f"{key}-{amounts}.in"
                create_file(filename, amounts, param["steps"])
                run_opalx(filename)

        for (key, param) in parameters.items():
            for amounts in param["amount"]:
                filename = f"{key}-{amounts}.stat"
                if not check_if_file_exists(f"opalx/{filename}"):
                    exit(0)
        MODE = "compare"

    if MODE == "local":
        for (key, param) in parameters.items():
            for amounts in param["amount"]:
                filename = f"{key}-{amounts}.in"
                create_file(filename, amounts, param["steps"])
                run_opalx(filename)
        MODE = "compare"
        
    if MODE == "compare":
        for (key, param) in parameters.items():
            for amounts in param["amount"]:
                filename = f"{key}-{amounts}.in"
                compare_data(filename, key)
        
        print("Plotting data...")
        calculate_statistics()
        plot_all()