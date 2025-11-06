import numpy as np
import scipy.constants as sc
from template import opalx_fodo_template
import matplotlib.pyplot as plt

FILE_OUT = "opalx/test.in"

N_PARTICLES   = 1e4 # in electrons amount
BUNCH_CHARGE  = sc.e * N_PARTICLES # C
INIT_ENERGY   = 1 # in GeV
DELTA_T       = 1e-10

LENGTH_DRIFT       = 2.5 # in m
LENGTH_QUADRUPOLE  = 0.5 # in m
FOCUSING_PARAMETER = 0.54102   # in 1/m^2 (alias k)
AMOUNT_OF_CELLS    = 100  

ENERGY_JOULES  = INIT_ENERGY * sc.e * 1e9 # in J
MOMENTUM = 1/sc.c * np.sqrt(ENERGY_JOULES**2 - (sc.m_e * sc.c**2)**2) # in kg m/s
BEAM_GRADIENT = FOCUSING_PARAMETER * (MOMENTUM / sc.e)  # T/m


def D(L): 
    return np.array([[1, L],
                     [0, 1]])

def Q(k,L):
    if k > 0: # for focusing quad
        j = np.sqrt(k)
        return np.array([[      np.cos(j * L), (1 / j) * np.sin(j * L)],
                         [- j * np.sin(j * L),           np.cos(j * L)]])
    else:  # for defocus
        j = np.sqrt(-k)
        return np.array([[    np.cosh(j * L), (1 / j) * np.sinh(j * L)],
                         [j * np.sinh(j * L),           np.cosh(j * L)]])

def round_sig(x, sig=2):
    if x == 0:
        return 0
    return np.round(x, sig - int(np.floor(np.log10(np.abs(x)))) - 1)

def propagate_twiss(alpha0, beta0, gamma0, s):
    M_total = np.eye(2)
    i = 0

    s_remaining = s
    while s_remaining > 0:
        if i == 0:
            L = min(s_remaining, LENGTH_QUADRUPOLE)
            M = Q(FOCUSING_PARAMETER, L)
        elif i == 1:
            L = min(s_remaining, LENGTH_DRIFT)
            M = D(L)
        elif i == 2:
            L = min(s_remaining, LENGTH_QUADRUPOLE)
            M = Q(-FOCUSING_PARAMETER, L)
        elif i == 3:
            L = min(s_remaining, LENGTH_DRIFT)
            M = D(L)
            i = -1
        s_remaining -= L
        i += 1
        M_total = M @ M_total  # pre-multiply

    M11, M12 = M_total[0,0], M_total[0,1]
    M21, M22 = M_total[1,0], M_total[1,1]

    beta  = M11**2 * beta0 - 2*M11*M12*alpha0 + M12**2 * gamma0
    alpha = -M11*M21*beta0 + (M11*M22 + M12*M21)*alpha0 - M12*M22*gamma0
    gamma = M21**2 * beta0 - 2*M21*M22*alpha0 + M22**2 * gamma0

    return np.array([alpha, beta, gamma])

def plot_theoretical_beta(alpha0, beta0, gamma0):
    x = np.linspace(0, (LENGTH_QUADRUPOLE + LENGTH_DRIFT) * 2 * AMOUNT_OF_CELLS, 1000)
    vals = np.array([propagate_twiss(alpha0, beta0, gamma0, s) for s in x])
    alpha, beta, gamma = vals[:, 0], vals[:, 1], vals[:, 2]

    plt.figure(figsize=[20, 10])
    plt.plot(x, alpha, label="alpha")
    plt.plot(x, beta,  label="beta")
    plt.plot(x, gamma,  label="gamma")
    
    plt.xticks(np.arange(AMOUNT_OF_CELLS) * (LENGTH_QUADRUPOLE + LENGTH_DRIFT) * 2)
    plt.legend()
    plt.grid()
    plt.savefig("theoretical vals")


M = Q(FOCUSING_PARAMETER, LENGTH_QUADRUPOLE) @ D(LENGTH_DRIFT) @ Q(-FOCUSING_PARAMETER, LENGTH_QUADRUPOLE) @ D(LENGTH_DRIFT) 
Tr = np.abs(np.trace(M))
phi = np.arccos(1/2 * Tr)
beta = M[0,1]/np.sin(phi)
alpha = (M[0,0] - np.cos(phi))/np.sin(phi)
gamma = (1 + alpha**2)/beta

alpha0 = alpha
beta0  = beta
gamma0 = gamma
if __name__ == "__main__":


    print("--- Theoretical values ---")
    if Tr < 2:
        print(f"The trace is |Tr| = {Tr} => stable")
    else:
        print(f"The trace is |Tr| = {Tr} => unstable")

    print(f"α = {alpha}")
    print(f"β = {beta} m")
    print(f"γ = {gamma} 1/m")
    print("")
    print(f"φ = {phi * 180/np.pi}°")
    print(f"G = {BEAM_GRADIENT} T/m")
    print("")
    plot_theoretical_beta(alpha, beta, gamma)
    print("Plotted theoretical values")

    LENGTH_OF_ALL_CELLS = AMOUNT_OF_CELLS * 2 * (LENGTH_DRIFT + LENGTH_QUADRUPOLE)
    velocity = MOMENTUM * sc.c**2 / ENERGY_JOULES
    time_needed = LENGTH_OF_ALL_CELLS / velocity  
    steps = np.ceil(time_needed / DELTA_T)

    min_distance = np.min([LENGTH_DRIFT, LENGTH_QUADRUPOLE])
    step_per_eleemnt = np.floor(min_distance / (velocity * DELTA_T))
    min_delta_t = min_distance/(velocity * 10)

    print("--- Simulation values ---")
    print(f"Length of all cells {LENGTH_OF_ALL_CELLS} m")
    print(f"Velocity {velocity/sc.c} c")
    if step_per_eleemnt < 10:
        print(f"Minmum step per element is {step_per_eleemnt:n} < 10. Δt was set to {round_sig(min_delta_t)}")
        DELTA_T = round_sig(min_delta_t)
        steps = np.ceil(time_needed / DELTA_T)
    print(f"Amount of steps needed {steps:n}")

    print("")

    with open(FILE_OUT, "w") as file:
        file.write(opalx_fodo_template(
            DELTA_T, BUNCH_CHARGE, N_PARTICLES, steps, INIT_ENERGY,  AMOUNT_OF_CELLS, LENGTH_QUADRUPOLE, LENGTH_DRIFT, BEAM_GRADIENT
        ))
        print(f"Written file to {FILE_OUT}")
