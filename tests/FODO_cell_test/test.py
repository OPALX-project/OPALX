import subprocess
import time

import sys
from opal_load_stat import load_dataset

import numpy as np
import pandas as pd
import os

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

opal = load_dataset(f"opalx/test.stat")
#opalx = load_dataset('./opalx', fname='multi.stat')

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

errors = {

}

def compare(stat1, stat2):
    amount = len(stat1.columns)
    passed = 0
    for c in stat2.columns:
        errors[c] = abs(stat1[c] - stat2[c]).to_numpy()
        if all(abs(stat1[c] - stat2[c]) < EPSILON) == False and c != "numParticles":
            print(f"\033[31mFailed {c} test: Max diff {max(abs(stat1[c] - stat2[c])):e}\033[0m")
            pass
        else:
            passed += 1

    ratio = passed/amount
    
    print(f"\033[35m --- Passed {ratio:.2%} of tests --- \033[0m")

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
    for (key, param) in errors.items():
        plt.figure(figsize=(10,5))

        plt.subplot(121)
        plt.semilogy(np.arange(1, len(param) + 1), param, label="Error")
        plt.xlabel("Steps")
        plt.ylabel("Error")
        plt.title(f"{key}")
        plt.legend()
        plt.grid()

        plt.subplot(122)
        plt.title(f"Data with $\Delta t =${opal['t'][0]:.0e}ns with {opal['numParticles'][0]} particles")
        plt.plot(np.arange(1, len(param) + 1), opalx[key], label="opalx")
        plt.plot(np.arange(1, len(param) + 1), opal[key], label="opal", ls= "--")
        plt.xlabel("Steps")
        plt.ylabel("Data")
        plt.legend()
        plt.grid()

        plt.savefig(f"{key}.png")
        plt.close()


#compare(opal, opalx)
#plot_all()

plt.figure(figsize=(7,5))
plt.title("RMS plot")
#plt.plot(opalx["s"], opalx["rms_x"], label="rms x opalx", color = "blue")
plt.plot(opal["s"], opal["rms_x"], ls = "--", label="rms x opal", color = "cyan")

#plt.plot(opalx["s"], opalx["rms_y"], label="rms y opalx", color ="red")
plt.plot(opal["s"], opal["rms_y"], ls = "--", label="rms y opal", color = "orange")

plt.grid()
plt.xlabel("s in [m]")
plt.ylabel("rms in [m]")
plt.legend()
plt.savefig("s-rms.png", dpi=300)
plt.close()


plt.figure(figsize=(7,5))
plt.title("Beta function")
plt.plot(opal["s"], opal["rms_x"]**2/opal["emit_x"])
plt.plot(opal["s"], opal["rms_y"]**2/opal["emit_y"])
plt.savefig("beta.png", dpi = 300)
plt.close()