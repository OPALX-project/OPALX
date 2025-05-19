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
OPAL_EXECUTABLE_FILE = "/data/user/binder_j/opal/build/src/opal"


# define the parameters
parameters = {
    "run1" : {
        "steps" : "30",
        "amount" : ["1e2", "1e3", "1e4"],
        "ref" : "ref.stat"
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
        print(f"Running {ORANGE}Opalx{WHITE} on dataset {BLUE + filename + WHITE}", end="", flush=True)

        start_time = time.time()
        result = subprocess.run([OPALX_EXECUTABLE_FILE, f"opalx/{filename}", "--info", "10"], capture_output=True, text=True, check=True)
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
        
def run_opal(filename):
    if os.path.isfile(f"opal/{filename.replace('in', 'stat')}"):
        print(f"{MAGENTA}Skipping execution with file {BLUE}{filename}{MAGENTA} for {ORANGE}Opal{WHITE}")
        return

    try:
        print(f"Running {ORANGE}Opal{WHITE} on dataset {BLUE + filename + WHITE}", end="", flush=True)
        start_time = time.time()
        result = subprocess.run([OPAL_EXECUTABLE_FILE, f"opal/{filename}", "--info", "10"], capture_output=True, text=True, check=True)
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
    print(f"Creating {BLUE+filename+WHITE}")
    with open(f"opal/{filename}", "w") as f:
        f.write(get_opal_string(amount, steps))
        f.close()
    
    with open(f"opalx/{filename}", "w") as f:
        f.write(get_opalx_string(amount, steps))
        f.close()

def compare(stat1, stat2, runname):
    amount = len(stat1.columns)
    passed = 0
    for c in stat2.columns:
        if c in parameters[runname]["data"].keys():
            # add the error if its observedd
            parameters[runname]["data"][c]["error"].append((abs(stat1[c] - stat2[c])).to_numpy())

        if all(abs(stat1[c] - stat2[c]) < EPSILON) == False:
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
    opal_filename  = filename.replace("in", "stat")
    print(f"Loading stats {BLUE}opalx/{opalx_filename + WHITE} and {BLUE}opal/{opal_filename + WHITE}")
    # load in the dataset
    opalxstat = load_dataset('opalx', fname=opalx_filename).dataframe
    opalstat = load_dataset('opal', fname=opal_filename).dataframe
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

for (key, param) in parameters.items():
    for amounts in param["amount"]:
        filename = f"{key}-{amounts}.in"
        create_file(filename, amounts, param["steps"])
        
        run_opal(filename)
        run_opalx(filename)

        compare_data(filename, key)
    pass

print("Plotting data...")
calculate_statistics()
plot_all()
print(f"{GREEN} --- Done --- {WHITE}")
#print(plotting_cols)
#plot_all()