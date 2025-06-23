import numpy as np
import scipy.constants as sc

BUNCH_CHARGE  = 1e-9 # C
N_PARTICLES   = 1e4 # in electrons amount
INIT_ENERGY   = 1e-3 # in GeV
BEAM_STRENGTH = .1 # T/m^2 

LENGTH_OF_QUADRUPOLE = 0.1 # in m
DRIFT_LENGTH = 0.5 # in m

DT = 1e-11

ELECTRON_MASS_GEV = 0.51099895000e-3  # GeV
ELECTRON_MASS_KG = sc.m_e  # kg

# Lorentz factors
gamma = (INIT_ENERGY + ELECTRON_MASS_GEV) / ELECTRON_MASS_GEV
beta = np.sqrt(1 - 1/gamma**2)

# Momentum in kg·m/s
p_si = gamma * beta * ELECTRON_MASS_KG * sc.c  

# Beam rigidity Bρ in Tesla·meters
beam_rigidity = p_si / sc.e  

# Focusing strength k [1/m^2]
k = BEAM_STRENGTH / beam_rigidity

# focal length f in [m]
focal_length = 1 / (k * LENGTH_OF_QUADRUPOLE) 

print(f"Focal length: {focal_length:.3f} m")

print(f"Maximum length of cell: {focal_length * 4:.3f} m")


POLE_LENGTH = LENGTH_OF_QUADRUPOLE # in m
amount_of_cells = 20

LENGTH_OF_ALL_CELLS = (POLE_LENGTH + DRIFT_LENGTH) * 2 * amount_of_cells # in m

# Velocity of the beam (beta * c)
velocity = beta * sc.c  

# Time to reach the distance
time_needed = LENGTH_OF_ALL_CELLS / velocity  

# Number of steps
steps = time_needed / DT

print(f"The length of all cells is {LENGTH_OF_ALL_CELLS:.3f} m\n")
print(f"Amount of steps needed: {np.ceil(steps):n}")

if (POLE_LENGTH + DRIFT_LENGTH) * 2 > focal_length * 4:
    print(f"Cell may be to big! Cell size: {(POLE_LENGTH + DRIFT_LENGTH) * 2} m")
    exit(0)
print(f"Cell size: {(POLE_LENGTH + DRIFT_LENGTH) * 2} m")



for i in range(amount_of_cells):
    start_of_cell = 2 * i * (POLE_LENGTH + DRIFT_LENGTH)
    print(f"FOCUS_POLE_{i}  : MULTIPOLE, L={POLE_LENGTH}, ELEMEDGE={start_of_cell:.3f}, KN={{0, {BEAM_STRENGTH}}};")
    print(f"FIRST_DRIFT_{i} : DRIFT    , L={DRIFT_LENGTH}, ELEMEDGE={start_of_cell + POLE_LENGTH:.3f};")
    print(f"DEFOCUS_POLE_{i}: MULTIPOLE, L={POLE_LENGTH}, ELEMEDGE={start_of_cell + POLE_LENGTH + DRIFT_LENGTH:.3f}, KN={{0, -{BEAM_STRENGTH}}};")
    print(f"SECOND_DRIFT_{i}: DRIFT    , L={DRIFT_LENGTH}, ELEMEDGE={start_of_cell + 2*POLE_LENGTH + DRIFT_LENGTH:.3f};")
    print("\n")

print("SOLine: Line = (")
for i in range(amount_of_cells - 1):
    print(f"FOCUS_POLE_{i}, FIRST_DRIFT_{i}, DEFOCUS_POLE_{i}, SECOND_DRIFT_{i},")
print(f"FOCUS_POLE_{amount_of_cells - 1}, FIRST_DRIFT_{amount_of_cells - 1}, DEFOCUS_POLE_{amount_of_cells - 1}, SECOND_DRIFT_{amount_of_cells - 1}")
print(");")