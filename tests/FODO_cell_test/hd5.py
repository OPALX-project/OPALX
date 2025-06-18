import sys
sys.path.append('/psi/home/adelmann/git/pyOPALTools')
from opal import load_dataset

import numpy as np
import pandas as pd
import h5py
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

px = []
x  = []

z = []

with h5py.File("opalx/multi.h5", "r") as f:
    for i in range(len(f.keys())):
        px.append((f[f"Step#{i}"]["px"][50]))
        x.append ((f[f"Step#{i}"]["x"] [50]))

        z.append ((f[f"Step#{i}"]["z"] [50]))

plt.plot(x, px, label="x-px")
plt.legend()
plt.grid()
plt.savefig("hd5out.png")
plt.close()

plt.figure(figsize=(20, 5))
steps = 1500
for i in range(5):
    plt.subplot(250 + i + 1)
    plt.plot(z[i * steps:(i + 1) * steps], x[i * steps:(i + 1) * steps] , label="z-x")
    plt.plot(z[i * steps:(i + 1) * steps], px[i * steps:(i + 1) * steps], label="z-px")
    plt.legend()
    plt.subplot(2, 5, i + 6)
    plt.plot(x[i * steps:(i + 1) * steps], px[i * steps:(i + 1) * steps], label="x-px")
    plt.legend()


plt.savefig("hd5out_allsteps.png", dpi=300)