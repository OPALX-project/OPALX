# import opal dataset
import sys
sys.path.append('/psi/home/adelmann/git/pyOPALTools')
from opal import load_dataset

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

EPSILON = 1e-9

#this is with the newly compiled version of opal
#opalstat = load_dataset('./reference/', fname='Solenoid-1.stat').dataframe

#opalxstat = load_dataset('/psi/home/adelmann/regression-tests/RegressionTests/Solenoid-1/', fname='Solenoid-1.stat').dataframe

opalxstat = load_dataset('/data/user/binder_j/opalx-elements/input-files/', fname='Solenoid-old.stat').dataframe
opalstat = load_dataset('./reference/', fname='Solenoid-0.stat').dataframe

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


def compare(stat1, stat2):
    amount = len(stat1.columns)
    passed = 0
    for c in stat1.columns:
        if all(abs(stat1[c] - stat2[c]) < EPSILON) == False:
            print(f"\033[31mFailed {c} test: Max diff {max(abs(stat1[c] - stat2[c])):.0e}\033[0m")  
        else:
            print(f"\033[34mPassed {c} test\033[0m")  
            passed += 1
    ratio = passed/amount
    
    print(f"\033[35m --- Passed {ratio:.2%} of tests --- \033[0m")
    

compare(opalstat, opalxstat)

plot_stat(opalstat)
plot_stat(opalxstat, title="opalx")
