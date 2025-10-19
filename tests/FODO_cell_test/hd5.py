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

def color(r, g, b, background = False):
    return f"\033[38;2;{r};{g};{b}m"

class H5Data:
    def __init__(self, x, y, z, px, py):
        self.x  = x
        self.y  = y
        self.z  = z
        self.px = px
        self.py = py

def load_h5(filename: str, particle_num = 1):
    print(f"Loading {filename}")
    with h5py.File(filename, "r") as f:
        steps = [f[f"Step#{i}"] for i in range(len(f.keys()))]
        px = np.array([s["px"][particle_num] for s in steps])
        py = np.array([s["py"][particle_num] for s in steps])
        x  = np.array([s["x"] [particle_num] for s in steps])
        y  = np.array([s["y"] [particle_num] for s in steps])
        z  = np.array([s["z"] [particle_num] for s in steps])

        print(f"Done loading...")
        return H5Data(x, y, z, px, py)

def plot_general(data: H5Data):
    plt.figure(figsize=(20, 5))
    plt.plot(REFERENCE_STAT["s"], data.x, label="x position", color="blue", ls="--")
    plt.plot(REFERENCE_STAT["s"], data.px, label="px momentum", color="blue")
    plt.plot(REFERENCE_STAT["s"], data.y, label="y position", color="red", ls="--")
    plt.plot(REFERENCE_STAT["s"], data.py, label="py momentum", color="red")
    plt.legend()
    plt.savefig("general_h5.png", bbox_inches="tight", dpi=300)

def fit_ellipse(x, px):
    # your points
    pts = np.column_stack((x, px))

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

    return (ellipse_x, ellipse_y)


def plot_ellipse(data: H5Data):
    # offsets to search for
    targets = np.arange(AMOUNT_OF_CELLS) * (DRIFT_LENGTH + LENGTH_OF_QUADRUPOLE) + (LENGTH_OF_QUADRUPOLE + DRIFT_LENGTH/2)

    # vectorized search
    indecis = np.searchsorted(REFERENCE_STAT["s"], targets, side='right')

    plt.figure(figsize=(12,10))
    plt.scatter(data.x[indecis], data.px[indecis], marker="x", label="x-component")
    plt.plot(*fit_ellipse(data.x[indecis], data.px[indecis]))
    plt.scatter(data.y[indecis], data.py[indecis], marker="x")
    plt.plot(*fit_ellipse(data.y[indecis], data.py[indecis]), label="y-component")
    plt.xlabel("Position space")
    plt.ylabel("Momentum space")
    plt.grid()
    plt.legend()
    plt.savefig("circle.png", bbox_inches="tight", dpi=300)
    plt.close()

def plot_all(data: H5Data):
    plot_general(data)
    plot_ellipse(data)

data = load_h5("opalx/test.h5")
plot_all(data)
