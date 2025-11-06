# import opal dataset
from opal_load_stat import load_dataset

import numpy as np
import os
import re

import matplotlib.pyplot as plt

from extract_time import parse_wall_avg

ROOT_FOLDER = "out"

FIT_PLOT = True

time_col = [
    "mainTimer",
    #"gatherInfoPartBunch",
    #"setSolver",
    #"samplingTime",
    #"GenParticles",
    "TIntegration1",
    "TIntegration2",
    #"External field eval",
    #"OrbThreader",
]

def compare_times():
    all_data = {}
    # get all times and put into a map
    for folder in os.listdir(ROOT_FOLDER):
        path = os.path.join(ROOT_FOLDER, folder)
        
        with open(os.path.join(path, "timing.dat"), "r") as file:
            lines = file.read()
            all_data[folder] = parse_wall_avg(lines.splitlines())
            print(f"Loaded {path}")
    
    reference_parameters = load_dataset(os.path.join(ROOT_FOLDER, os.listdir(ROOT_FOLDER)[0], "template.stat"))

    check = re.compile(r"\D*(\d+)\D*(\d+)")

    print("---")
    plt.figure(figsize=(10, 8))
    for col in time_col:
        data = []
        for name in all_data.keys():
            # extract the amount of nodes and threads that were used for a run
            nodes   = int(check.match(name)[1])
            threads = int(check.match(name)[2])
            data.append([nodes * threads, all_data[name][col]])

        # sort them
        data = np.array(sorted(data, key=lambda x: x[0]))
        
        x = data[:,0]
        y = data[:,1]

        print(f"Plotting {col}")
        ((m, b), cov) = np.polyfit(np.log10(x), np.log10(y), 1, cov = True)
        m_err, b_err = np.sqrt(np.diag(cov))

        # and plot
        fitted = 10**b * x**m  # convert back from log-log

        # calculate the envelope
        x_log = np.log10(x)
        z_fit = m * x_log + b
        y_fit = 10**z_fit

        var_m  = cov[0, 0]
        var_b  = cov[1, 1]
        cov_mb = cov[0, 1]
        var_z = var_b + (x_log**2) * var_m + 2 * x_log * cov_mb
        sigma_z = np.sqrt(var_z)

        y_lo = 10**(z_fit - sigma_z)
        y_hi = 10**(z_fit + sigma_z)

        
        line, = plt.plot(x, y, ls="--", marker="D", label="" if FIT_PLOT else col, zorder=100)
        color = line.get_color()
        if FIT_PLOT:
            plt.plot(x, fitted, ls="-", label=f"{col} (slope = {m:.2f})", alpha=0.5, color=color)
            plt.fill_between(x, y_lo, y_hi, color=color, alpha=0.1)
    
    # using the reference solution print out more data

    plt.yscale("log", base=10)   # log base 10 on y-axis
    plt.xscale("log", base=2)    # log base 2 on x-axis
    plt.xlabel("Amount of OpenMP threads")
    plt.ylabel("Time in [s]")
    plt.legend()
    plt.grid()
    plt.title(f"Speedup with {reference_parameters['numParticles'][0]:.0e} particles, {len(reference_parameters['numParticles'])} steps @ $\Delta t$ = {reference_parameters['t'][0]} ns")
    plt.savefig("times.png", bbox_inches="tight", dpi=300)

if __name__ == "__main__":
    compare_times()
    print("Done")
