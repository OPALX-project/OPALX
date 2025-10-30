import sys
import scipy.constants as sc
from opal_load_stat import load_dataset

import generate #used to get the parameters

OPALX_STAT = load_dataset("opalx/test.stat")
OPAL_STAT = load_dataset("opal/test.stat")

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import h5py
import os

import matplotlib.pyplot as plt

def color(r, g, b, background = False):
    return f"\033[38;2;{r};{g};{b}m"


def plot_general():
    # Assuming data.x, data.px, data.y, and data.py are your data arrays
    fig, ax1 = plt.subplots(figsize=(20, 5))

    # Plot position data on the first y-axis
    ax1.plot(OPALX_STAT["s"], OPALX_STAT["ref_y"], label="OPALX x position", color="blue", ls="--")
    ax1.plot(OPALX_STAT["s"], OPALX_STAT["ref_y"], label="OPALX y position", color="red", ls="--")
    ax1.set_xlabel("z position in the beam line [m]")
    ax1.set_ylabel("Position [m]")
    ax1.tick_params(axis='y')
    ax1.grid()

    # Create a second y-axis for momentum data
    ax2 = ax1.twinx()
    ax2.plot(OPALX_STAT["s"], OPALX_STAT["ref_px"], label="px momentum", color="blue")
    ax2.plot(OPALX_STAT["s"], OPALX_STAT["ref_py"], label="py momentum", color="red")
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

    # covariance matrix
    S = np.cov(np.vstack([x, y]))
    sig_xx, sig_xpx, sig_pxpx = S[0,0], S[0,1], S[1,1]

    # geometric emittance
    eps = np.sqrt(sig_xx * sig_pxpx - sig_xpx**2)

    # Twiss parameters
    beta  = sig_xx / eps
    alpha = -sig_xpx / eps
    gamma = (1 + alpha**2) / beta

    print(alpha)
    print(beta)
    print(gamma)
    # draw ellipse
    phi = np.linspace(0, 2*np.pi, 300)
    ellipse_x = a*np.cos(phi)*ct - b*np.sin(phi)*st
    ellipse_y = a*np.cos(phi)*st + b*np.sin(phi)*ct

    return ellipse_x, ellipse_y


def plot_ellipse():
    # offsets to search for
    targets = np.arange(generate.AMOUNT_OF_CELLS) * 2 * (generate.LENGTH_DRIFT + generate.LENGTH_QUADRUPOLE) + (generate.LENGTH_QUADRUPOLE + generate.LENGTH_DRIFT / 2)

    # vectorized search
    indecis = np.searchsorted(OPALX_STAT["s"], targets, side='right')
    plt.figure(figsize=(12,10))

    with h5py.File("opalx/test.h5", "r") as f:
        steps = [f[f"Step#{i}"] for i in indecis]
        px = np.array([s["px"][1000] for s in steps])
        py = np.array([s["py"][1000] for s in steps])
        x  = np.array([s["x"] [1000] for s in steps])
        y  = np.array([s["y"] [1000] for s in steps])
    plt.scatter(x, px, marker="x", label="OPAL-X x-px-component")
    plt.plot(*fit_ellipse(x, px), label="Fitted ellipse for the x component")
    plt.scatter(y, py, marker="x", label="OPAL-X y-py-component")
    plt.plot(*fit_ellipse(y, py), label="Fitted ellipse for the x component")

    plt.xlabel("Position space [m]")
    plt.ylabel("Momentum space")
    plt.title("OPAL-X Particle position-momentum plot at even distances")
    plt.grid()
    plt.legend()
    plt.savefig("circle.png", bbox_inches="tight", dpi=300)
    plt.close()

def plot_beta():
    plt.figure(figsize=(10,10))
    
    fodo = np.arange(generate.AMOUNT_OF_CELLS) * (generate.LENGTH_DRIFT + generate.LENGTH_QUADRUPOLE) * 2
    
    def calc_beta(data):
        s =      np.asarray(data["s"])
        rms_x =  np.asarray(data["rms_x"])
        rms_y =  np.asarray(data["rms_y"])
        # normalized emittance
        emit_x = np.asarray(data["emit_x"])
        emit_y = np.asarray(data["emit_y"])

        energy_joules = np.asarray(data["energy"]) * sc.e * 1e6
        gamma = energy_joules/(sc.c**2 * sc.m_e)
        beta = np.sqrt(1 - 1/gamma**2)

        emit_x_geo = emit_x / (gamma * beta) + 1e-9
        emit_y_geo = emit_y / (gamma * beta) + 1e-9
        # beta functions
        beta_x = rms_x**2 / emit_x_geo
        beta_y = rms_y**2 / emit_y_geo

        return (beta_x, beta_y)
        

    #theorey
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 14), sharex=True)
    # --- Beta function plot ---
    beta_x, beta_y = calc_beta(OPAL_STAT)
    opalx_beta_x, opalx_beta_y = calc_beta(OPALX_STAT)

    s = OPAL_STAT["s"]
    vals = np.array([generate.propagate_twiss(generate.alpha0, generate.beta0, generate.gamma0, x) for x in s])
    _, beta, _ = vals[:, 0], vals[:, 1], vals[:, 2]
    ax1.plot(s, beta_x, label="OPAL βx", ls="--", color="red")
    ax1.plot(s, beta_y, label="OPAL βy", ls="--", color="blue")
    ax1.plot(OPALX_STAT["s"], opalx_beta_x, label="OPAL-X βx", color="red")
    ax1.plot(OPALX_STAT["s"], opalx_beta_y, label="OPAL-X βy", color="blue")
    ax1.plot(s, beta, label="Theoretical β", ls="dashdot", color="orange")

    ax1.set_ylabel("β [m]")
    ax1.set_title("Beta Function along FODO Cell")
    ax1.set_ylim([0,120])
    ax1.grid(True)
    ax1.legend()
    
    # offsets to search for
    targets = np.arange(generate.AMOUNT_OF_CELLS) * 2 * (generate.LENGTH_DRIFT + generate.LENGTH_QUADRUPOLE) + (generate.LENGTH_QUADRUPOLE + generate.LENGTH_DRIFT / 2)

    # vectorized search
    indecis = np.searchsorted(OPALX_STAT["s"], targets, side='right')

    # --- Error subplot (difference between OPAL-X and OPAL) ---
    err_x = np.abs(opalx_beta_x[1:] - beta_x) / beta_x
    err_y = np.abs(opalx_beta_y[1:] - beta_y) / beta_y

    ax2.plot(s, err_x * 100, label="βx Error", color="red")
    ax2.plot(s, err_y * 100, label="βy Error", color="blue")
    ax2.set_ylim([0,100])
    ax2.set_xlabel("s [m]")
    ax2.set_ylabel("Relative error Δβ/β [%]")
    ax2.grid(True)
    ax2.legend()

    # plot
    plt.xlim([300,600])
    plt.grid(True)
    plt.legend()
    plt.savefig("beta.png", bbox_inches="tight", dpi=300)

def plot_all():
    #plot_ellipse()
    #plot_general()
    plot_beta()

plot_all()
