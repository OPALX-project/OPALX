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
    x =  np.asarray(REFERENCE_STAT["ref_x"])
    y =  np.asarray(REFERENCE_STAT["ref_y"])
    z =  np.asarray(REFERENCE_STAT["ref_z"])
    px = np.asarray(REFERENCE_STAT["ref_px"])
    py = np.asarray(REFERENCE_STAT["ref_py"])
    return H5Data(x, y, z, px, py)

def plot_general(data: H5Data):
    # Assuming data.x, data.px, data.y, and data.py are your data arrays
    fig, ax1 = plt.subplots(figsize=(20, 5))

    # Plot position data on the first y-axis
    ax1.plot(REFERENCE_STAT["s"], data.x, label="x position", color="blue", ls="--")
    ax1.plot(REFERENCE_STAT["s"], data.y, label="y position", color="red", ls="--")
    ax1.set_xlabel("z position in the beam line [m]")
    ax1.set_ylabel("Position [m]")
    ax1.tick_params(axis='y')
    ax1.grid()

    # Create a second y-axis for momentum data
    ax2 = ax1.twinx()
    ax2.plot(REFERENCE_STAT["s"], data.px, label="px momentum", color="blue")
    ax2.plot(REFERENCE_STAT["s"], data.py, label="py momentum", color="red")
    ax2.set_ylabel("Momentum")
    ax2.tick_params(axis='y')
    ax2.grid(ls =":")

    # Add legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="best")
    ax1.set_title("Reference particle position and momentum vs distance")

    plt.savefig("general_h5.png", bbox_inches="tight", dpi=300)

def fit_ellipse(x, y):
    # build design matrix for Ax^2 + Bxy + Cy^2 + 1 = 0
    D = np.vstack([x**2, x*y, y**2]).T
    d = -np.ones_like(x)

    # solve least squares
    A, B, C = np.linalg.lstsq(D, d, rcond=None)[0]

    # extract parameters
    theta = 0.5 * np.arctan2(B, A - C)
    ct, st = np.cos(theta), np.sin(theta)
    Ap = A*ct**2 + B*ct*st + C*st**2
    Cp = A*st**2 - B*ct*st + C*ct**2
    a = np.sqrt(-1 / Ap)
    b = np.sqrt(-1 / Cp)

    # draw ellipse
    phi = np.linspace(0, 2*np.pi, 300)
    ellipse_x = a*np.cos(phi)*ct - b*np.sin(phi)*st
    ellipse_y = a*np.cos(phi)*st + b*np.sin(phi)*ct

    return ellipse_x, ellipse_y


def plot_ellipse(data: H5Data):
    # offsets to search for
    targets = np.arange(AMOUNT_OF_CELLS) * (DRIFT_LENGTH + LENGTH_OF_QUADRUPOLE) + (LENGTH_OF_QUADRUPOLE + DRIFT_LENGTH/2)

    # vectorized search
    indecis = np.searchsorted(REFERENCE_STAT["s"], targets, side='right')

    plt.figure(figsize=(12,10))
    plt.scatter(data.x[indecis], data.px[indecis], marker="x", label="x-px-component")
    plt.plot(*fit_ellipse(data.x[indecis], data.px[indecis]), label="Fitted ellipse for the x component")
    plt.scatter(data.y[indecis], data.py[indecis], marker="x", label="y-py-component")
    plt.plot(*fit_ellipse(data.y[indecis], data.py[indecis]), label="Fitted ellipse for the y component")
    plt.xlabel("Position space [m]")
    plt.ylabel("Momentum space")
    plt.title("Reference particle position-momentum plot at even distances")
    plt.grid()
    plt.legend()
    plt.savefig("circle.png", bbox_inches="tight", dpi=300)
    plt.close()

def plot_beta(data: H5Data):
    plt.figure(figsize=(10,10))
    
    s =      REFERENCE_STAT["s"]
    rms_x =  REFERENCE_STAT["rms_x"]
    rms_y =  REFERENCE_STAT["rms_y"]
    emit_x = REFERENCE_STAT["emit_x"]
    emit_y = REFERENCE_STAT["emit_y"]

    # beta functions
    beta_x = rms_x**2 / emit_x
    beta_y = rms_y**2 / emit_y

    # plot
    #plt.plot(s, beta_x, label="βx")
    #plt.plot(s, beta_y, label="βy")
    plt.plot(s, emit_x, label="emi_x")
    plt.xlabel("s [m]")
    plt.ylabel("β [m]")
    plt.title("Beta Function along FODO Cell")
    plt.grid(True)
    plt.legend()
    plt.savefig("beta.png", bbox_inches="tight", dpi=300)

def plot_all(data: H5Data):
    plot_general(data)
    plot_ellipse(data)
    plot_beta(data)

data = load_h5("opalx/test.h5")
plot_all(data)
