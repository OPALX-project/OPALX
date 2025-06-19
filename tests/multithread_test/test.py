# import opal dataset
import sys
sys.path.append('/psi/home/adelmann/git/pyOPALTools')
from opal import load_dataset

import numpy as np
import os
import re
import pandas as pd

import matplotlib.pyplot as plt

from extract_time import parse_wall_avg

EPSILON = 1e-9

ROOT_FOLDER = "out"
all_folders = os.listdir(ROOT_FOLDER)
REFERENCE_FOLDER = "n1t1" # single thread should be the "ref" solution
all_folders = list(filter(lambda x : x != REFERENCE_FOLDER, all_folders))

REFERENCE_SOLUTION = load_dataset(os.path.join(ROOT_FOLDER, REFERENCE_FOLDER), fname="template.stat").dataframe

failed_col = []

time_col = [
    "mainTimer",
    "gatherInfoPartBunch",
    "setSolver",
    "samplingTime",
    "GenParticles",
    "TIntegration1",
    "TIntegration2",
    "External field eval",
    "OrbThreader",
]

def compare_all_cols():
    all_data = {}
    for folder_name in all_folders:
        print(f" --- {folder_name} ---")
        data = load_dataset(os.path.join(ROOT_FOLDER, folder_name), fname="template.stat").dataframe
        all_data[folder_name] = data

        # iterate over all coloumns
        compare(REFERENCE_SOLUTION, data)

    # iterate over all cols
    
    for c in REFERENCE_SOLUTION.columns:
        if not (c in failed_col):
            continue

        plt.figure(figsize=(10, 5))
        plt.subplot(121)
        x = np.arange(1, len(REFERENCE_SOLUTION[c]) + 1)

        # now go over multithreaded solution
        for d in all_data.keys():
            plt.plot(x, all_data[d][c], label=f"{d}")

        plt.plot(x, REFERENCE_SOLUTION[c], label="Reference", ls="--")
        plt.legend()
        plt.xlabel("Steps")
        plt.title("Raw data")
        plt.grid()

        plt.subplot(122)
        for d in all_data.keys():
            plt.semilogy(x, np.abs(all_data[d][c] - REFERENCE_SOLUTION[c]), label=f"{d}")
        plt.xlabel("Steps")
        plt.title("Error")
        plt.grid()
        

        plt.savefig(f"{c}.png")
        plt.close()

def compare_times():
    all_data = {}
    # get all times and put into a map
    for folder in os.listdir(ROOT_FOLDER):
        path = os.path.join(ROOT_FOLDER, folder)
        
        with open(os.path.join(path, "timing.dat"), "r") as file:
            lines = file.read()
            all_data[folder] = parse_wall_avg(lines.splitlines())
    

    check = re.compile(r"\D*(\d+)\D*(\d+)")

    plt.figure(figsize=(12, 10))
    for col in time_col:
        data = []
        for name in all_data.keys():
            # extract the amount of nodes and threads that were used for a run
            nodes   = int(check.match(name)[1])
            threads = int(check.match(name)[2])
            data.append([nodes * threads, all_data[name][col]])

        # sort them
        data = np.array(sorted(data, key=lambda x: x[0]))

        # and plot
        plt.plot(data[:,0], data[:,1], ls="--", marker="D", label=col, zorder=100)
    
    # using the reference solution print out more data

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel("Amount of threads")
    plt.ylabel("Time in [s]")
    plt.legend()
    plt.grid()
    plt.title(f"Speedup with {REFERENCE_SOLUTION['numParticles'][0]:.0e} particles, {len(REFERENCE_SOLUTION['numParticles'])} steps @ $\Delta t$ = {REFERENCE_SOLUTION['t'][0]} ns")
    plt.savefig("times.png")

        
        

def compare(stat1, stat2):
    amount = len(stat1.columns)
    passed = 0

    it = 0
    for c in stat1.columns:
        it += 1
        if all(abs(stat1[c] - stat2[c]) < EPSILON) == False:
            failed_col.append(c)
            print(f"\033[31mFailed {c} test: Max diff {max(abs(stat1[c] - stat2[c])):.0e}\033[0m")  
        else:
            #print(f"\033[34mPassed {c} test\033[0m")  
            passed += 1
    ratio = passed/amount
    
    print(f"\033[35m --- Passed {ratio:.2%} of tests --- \033[0m")

# not needed anymore   
#compare_all_cols()
compare_times()