import sys
from opal_load_stat import load_dataset

REFERENCE_STAT = load_dataset("opalx/test.stat")
LENGTH_OF_QUADRUPOLE = 0.005 # in m
DRIFT_LENGTH = 0.05 # in m
AMOUNT_OF_CELLS = 50

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import h5py
import os

import matplotlib.pyplot as plt

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
py = []
x  = []
y  = []
z = []
indecis = []

print()

with h5py.File("opalx/test.h5", "r") as f:
    print(len(f.keys()))
    for i in range(len(f.keys())):
        px.append(f[f"Step#{i}"]["px"][1500])
        py.append(f[f"Step#{i}"]["py"][1500])
        x.append (f[f"Step#{i}"]["x"] [1500])
        y.append (f[f"Step#{i}"]["y"] [1500])
        z.append (f[f"Step#{i}"]["z"] [1500])

x  = np.asarray(x)
px = np.asarray(px)
y  = np.asarray(y)
py = np.asarray(py)

# offsets to search for
targets = np.arange(AMOUNT_OF_CELLS) * (DRIFT_LENGTH + LENGTH_OF_QUADRUPOLE) + (LENGTH_OF_QUADRUPOLE + DRIFT_LENGTH/2)

# vectorized search
indecis = np.searchsorted(REFERENCE_STAT["s"], targets, side='right')

print(indecis)
plt.figure(figsize=(12, 10))

plt.scatter(x[indecis], px[indecis])
plt.scatter(y[indecis], py[indecis])
plt.savefig("circle.png")
plt.close()

# your points
pts = np.column_stack((x[indecis], px[indecis]))

# algebraic ellipse fit
def ellipse_res(params, x, y):
    xc, yc, a, b, theta = params
    ct, st = np.cos(theta), np.sin(theta)
    xr = ct*(x-xc) + st*(y-yc)
    yr = -st*(x-xc) + ct*(y-yc)
    return (xr/a)**2 + (yr/b)**2 - 1

x0 = np.mean(pts[:,0])
y0 = np.mean(pts[:,1])
a0 = (pts[:,0].max() - pts[:,0].min())/2
b0 = (pts[:,1].max() - pts[:,1].min())/2
theta0 = 0

res = least_squares(ellipse_res, x0=[x0,y0,a0,b0,theta0],
                    args=(pts[:,0], pts[:,1]))
xc, yc, a, b, theta = res.x

# draw ellipse
phi = np.linspace(0, 2*np.pi, 300)
ct, st = np.cos(theta), np.sin(theta)
ellipse_x = xc + a*np.cos(phi)*ct - b*np.sin(phi)*st
ellipse_y = yc + a*np.cos(phi)*st + b*np.sin(phi)*ct

plt.figure(figsize=(12,10))
plt.scatter(pts[:,0], pts[:,1])
plt.plot(ellipse_x, ellipse_y)
plt.savefig("circle.png")
plt.close()

plt.figure(figsize=(12, 10))
plt.plot(x, px, label="x-px")
plt.plot(y, py, label="y-py")
plt.legend()
plt.grid()
plt.savefig("hd5out.png")
plt.close()

plt.figure(figsize=(20, 5))
plt.plot(REFERENCE_STAT["s"], x)
plt.plot(REFERENCE_STAT["s"], px)
plt.plot(REFERENCE_STAT["s"], y)
plt.plot(REFERENCE_STAT["s"], py)
plt.savefig("hd5out2.png")
# steps = 1500
# for i in range(5):
#     plt.subplot(250 + i + 1)
#     plt.plot(z[i * steps:(i + 1) * steps], x[i * steps:(i + 1) * steps] , label="z-x")
#     plt.plot(z[i * steps:(i + 1) * steps], px[i * steps:(i + 1) * steps], label="z-px")
#     plt.legend()
#     plt.subplot(2, 5, i + 6)
#     plt.plot(x[i * steps:(i + 1) * steps], px[i * steps:(i + 1) * steps], label="x-px")
#     plt.legend()


# plt.savefig("hd5out_allsteps.png", dpi=300)