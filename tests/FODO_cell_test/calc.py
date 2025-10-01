import numpy as np
import scipy.constants as sc

N_PARTICLES   = 1e4 # in electrons amount
BUNCH_CHARGE  = sc.e * N_PARTICLES # C
INIT_ENERGY   = 1e-3 # in GeV


BEAM_STRENGTH = 0.00257 # T/m 
DT = 2e-11
LENGTH_OF_QUADRUPOLE = 0.25 # in m
DRIFT_LENGTH = 2.5 # in m
amount_of_cells = 9



ELECTRON_MASS_GEV = 0.51099895000e-3  # GeV
ELECTRON_MASS_KG = sc.m_e  # kg

# Lorentz factors
gamma = (INIT_ENERGY + ELECTRON_MASS_GEV) / ELECTRON_MASS_GEV
beta = np.sqrt(1 - 1/gamma**2)

# Momentum in kg·m/s
p_si = gamma * beta * ELECTRON_MASS_KG * sc.c  

# Beam rigidity Bρ in Tesla·meters
beam_rigidity = p_si / (BUNCH_CHARGE / N_PARTICLES)

# Focusing strength k [1/m^2]
k = BEAM_STRENGTH / beam_rigidity

# focal length f in [m]
focal_length = 1 / (k * LENGTH_OF_QUADRUPOLE) 

print(f"Focal length: {focal_length:.3f} m")
print(f"Beta: {beta} c")
print(f"K: {k:4f} 1/m²")

print(f"Maximum length of cell: {focal_length * 4:.3f} m")


POLE_LENGTH = LENGTH_OF_QUADRUPOLE # in m

LENGTH_OF_ALL_CELLS = (POLE_LENGTH + DRIFT_LENGTH) * 2 * amount_of_cells # in m

# Velocity of the beam (beta * c)
velocity = beta * sc.c  

# Time to reach the distance
time_needed = LENGTH_OF_ALL_CELLS / velocity  

# Number of steps
steps = np.ceil(time_needed / DT)

print(f"The length of all cells is {LENGTH_OF_ALL_CELLS:.3f} m\n")
print(f"Amount of steps needed: {np.ceil(steps):n}")

if (POLE_LENGTH + DRIFT_LENGTH) * 2 > focal_length * 4:
    print(f"Cell may be to big! Cell size: {(POLE_LENGTH + DRIFT_LENGTH) * 2} m")
    #exit(-1)
print(f"Cell size: {(POLE_LENGTH + DRIFT_LENGTH) * 2} m")

def opal_template(s):
    global DT
    global BUNCH_CHARGE
    global N_PARTICLES
    global steps
    return f"""OPTION, PSDUMPFREQ = 1;   // 6d data written every 300000 time steps (h5).
OPTION, STATDUMPFREQ = 1;  // Beam Stats written every 10 time steps (stat).
OPTION, BOUNDPDESTROYFQ=10; // Delete lost particles, if any
OPTION, AUTOPHASE=4;        // Autophase is on, and phase of max energy
                            // gain will be found automatically for cavities
Option, VERSION=10900;

Title, string="Solenoid-1 test";

//----------------------------------------------------------------------------
//Global Parameters

REAL rf_freq             = 1.3e3;     //RF frequency. (Hz)
REAL n_particles         = {N_PARTICLES:n};      //Number of particles in simulation.
REAL beam_bunch_charge   = {BUNCH_CHARGE};      //Charge of bunch. (C)

//Initial Momentum Calculation
REAL Edes    = 1e-3;  // initial energy in GeV
REAL gamma   = (Edes+EMASS)/EMASS;
REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;    //inital z momentum

//Printing initial energy and momentum to terminal output.
value , {{Edes, P0, OPALVERSION}};

{s}

Dist: DISTRIBUTION, TYPE = GAUSS,
        SIGMAX = 0.001,
        SIGMAY = 0.001,
	SIGMAZ = 0.001;

//----------------------------------------------------------------------------
// Define Field solvers
// The mesh sizes should be a factor of 2
// for most efficient space charge calculation.

FS_SC: Fieldsolver, FSTYPE = NONE,
            MX = 16, MY = 16, MT = 16, // SC grid size is 32^3
            PARFFTX = false,
            PARFFTY = false,
            PARFFTT = true,  // parallel in the z direction only
            BCFFTX = open,
            BCFFTY = open,
            BCFFTT = open,
            BBOXINCR = 1,
            GREENSF = INTEGRATED;

//----------------------------------------------------------------------------
// Electron Beam Definition

BEAM1:  BEAM, PARTICLE = ELECTRON, pc = P0, NPART = n_particles,
        BFREQ = rf_freq, BCURRENT = beam_bunch_charge * rf_freq * 1e6, CHARGE = -1;

//----------------------------------------------------------------------------
// Simulate the beamline using TRACK and RUN.
// Note, different time steps are set based on the z location in the beam line.
// In the case below, 1.0e-13 is used for 0.0 to 0.4 m,
// and 3.0e-12 is used from 0.4 to 5 m.

TRACK, LINE = SOLine, BEAM = BEAM1, MAXSTEPS = {steps:n}, DT = {DT};
RUN, METHOD = "PARALLEL-T", BEAM = BEAM1,
    FIELDSOLVER = FS_SC, DISTRIBUTION = Dist;
ENDTRACK;
Quit;
"""

out_string = ""

for i in range(amount_of_cells):
    start_of_cell = 2 * i * (POLE_LENGTH + DRIFT_LENGTH)
    out_string += (f"FOCUS_POLE_{i}  : MULTIPOLE, L={POLE_LENGTH}, ELEMEDGE={start_of_cell:.3f}, KN={{0, {BEAM_STRENGTH}}};\n")
    out_string += (f"FIRST_DRIFT_{i} : DRIFT    , L={DRIFT_LENGTH}, ELEMEDGE={start_of_cell + POLE_LENGTH:.3f};\n")
    out_string += (f"DEFOCUS_POLE_{i}: MULTIPOLE, L={POLE_LENGTH}, ELEMEDGE={start_of_cell + POLE_LENGTH + DRIFT_LENGTH:.3f}, KN={{0, -{BEAM_STRENGTH}}};\n")
    out_string += (f"SECOND_DRIFT_{i}: DRIFT    , L={DRIFT_LENGTH}, ELEMEDGE={start_of_cell + 2*POLE_LENGTH + DRIFT_LENGTH:.3f};\n")
    out_string += ("\n")

out_string += ("SOLine: Line = (")
for i in range(amount_of_cells - 1):
    out_string += (f"FOCUS_POLE_{i}, FIRST_DRIFT_{i}, DEFOCUS_POLE_{i}, SECOND_DRIFT_{i},\n")
out_string += (f"FOCUS_POLE_{amount_of_cells - 1}, FIRST_DRIFT_{amount_of_cells - 1}, DEFOCUS_POLE_{amount_of_cells - 1}, SECOND_DRIFT_{amount_of_cells - 1}\n")
out_string += (");")

with open("opal/test.in", "w") as file:
    file.write(opal_template(out_string))

print("Written file to test.in")