# import opal dataset
import sys
sys.path.append('/psi/home/adelmann/git/pyOPALTools')
from opal import load_dataset

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

EPSILON = 1e-9

errors = {}

multi_thread = load_dataset('./single', fname='single.stat').dataframe
single_thread = load_dataset('./multi', fname='multi.stat').dataframe

def plot_stat(stat, title="opal"):
    plt.figure(figsize=(10, 10))
    plt.title(title)
    plt.subplot(2,2,1)
    #plt.ylim((-0.125, 0.125))
    
    plt.title("rms")
    plt.xlabel('s [m]')
    plt.ylabel('rms x [m]')
    plt.grid()
    plt.plot(stat['s'],100*stat['rms_y'],label=f"y *100 {title}")
    plt.plot(stat['s'],100*stat['rms_x'],label=f"x *100 {title}")
    plt.plot(stat['s'],100*stat['rms_s'],label=f"z *100 {title}")
    plt.legend()
    
    plt.subplot(2, 2, 2)
    plt.grid()
    plt.plot(stat['s'],stat['Bz_ref'] * 100,label=f"Bz {title}")
    plt.xlabel('s [m]')
    plt.ylabel('B [T]')
    
    plt.subplot(2, 2, 3)
    plt.grid()
    plt.title("rms_p")
    plt.plot(stat['s'],100*stat['rms_py'],label=f"y *100 {title}")
    plt.plot(stat['s'],100*stat['rms_px'],label=f"x *100 {title}")
    plt.plot(stat['s'],100*stat['rms_ps'],label=f"z *100 {title}")
    
    plt.subplot(2, 2, 4)
    plt.grid()
    plt.title("emit")
    plt.plot(stat['s'],100*stat['emit_x'],label=f"y *100 {title}")
    plt.plot(stat['s'],100*stat['emit_y'],label=f"x *100 {title}")
    plt.plot(stat['s'],100*stat['emit_s'],label=f"z *100 {title}")

    plt.savefig(f"{title}.png", dpi=300)
    print("Plotted all stats")

def plot_errors():
    first_y = errors["emit_x"]
    secon_y = errors["rms_x"]
    x = np.linspace(0, len(first_y) - 1, len(first_y))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

    ax1.set_title("500 Steps at $\\Delta t = 10^{-11}$")
    ax1_twin = ax1.twinx()
    ax2_twin = ax2.twinx()

    line1,  = ax1.plot(x, multi_thread["rms_x"], label="rms_x", color = "blue")
    line1_ , = ax1.plot(x, single_thread["rms_x"], linestyle = ":", color = "green")
    ax1.tick_params(axis="y", color="blue")
    ax1.set_ylabel("rms_x", color="blue")

    line2, = ax1_twin.plot(x, multi_thread["emit_x"], label = "emit_x", color="red")
    line2_, = ax1_twin.plot(x, single_thread["emit_x"], linestyle = ":", color = "orange")

    ax1_twin.tick_params(axis="y", color="red")
    ax1_twin.set_ylabel('Emit_x', color='red', rotation =270)
    ax1.grid()

    ax1.legend([line1, line1_, line2, line2_], ["rms_x multiple threads", "rms_x single thread", "emit_x multiple threads", "emit_x single thread"], loc=9) 
    ax1.grid(color="blue", alpha=0.2)
    ax1_twin.grid(color="red", alpha=0.2)

    line1, = ax2.semilogy(x, secon_y, color = "blue", label="rms_x error")
    line2, = ax2_twin.semilogy(x, first_y, color="red", label="emit_x error")
    ax2.tick_params(axis="y", color="blue")
    ax2.set_ylabel("rms_x", color="blue")
    ax2_twin.tick_params(axis="y", color="red")
    ax2_twin.set_ylabel('Emit_x', color='red', rotation =270)
    ax2.legend([line1, line2], ["rms_x error", "emit_x error"], loc=7) 

    ax2.grid(color="blue", alpha=0.2)
    ax2_twin.grid(color="red", alpha=0.2)

    
    ax1.set_xlabel("Steps")

    plt.savefig("error.png", bbox_inches = "tight")


def compare(stat1, stat2):
    amount = len(stat1.columns)
    passed = 0

    it = 0
    for c in stat1.columns:
        it += 1
        errors[c] = abs(stat1[c] - stat2[c]).to_numpy()
        if all(abs(stat1[c] - stat2[c]) < EPSILON) == False:
            print(f"\033[31mFailed {c} test: Max diff {max(abs(stat1[c] - stat2[c])):.0e}\033[0m")  
        else:
            print(f"\033[34mPassed {c} test\033[0m")  
            passed += 1
    ratio = passed/amount
    
    print(f"\033[35m --- Passed {ratio:.2%} of tests --- \033[0m")
    

compare(single_thread, multi_thread)
plot_errors()