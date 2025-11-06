def opalx_fodo_template(dt, bunch_charge, n_parts, steps, energy,
                        amount_of_cells, POLE_LENGTH, DRIFT_LENGTH, BEAM_STRENGTH):
    input_file = lambda s : f"""OPTION, PSDUMPFREQ = 1;   // 6d data written every 300000 time steps (h5).
OPTION, STATDUMPFREQ = 1;  // Beam Stats written every 10 time steps (stat).
OPTION, BOUNDPDESTROYFQ=10; // Delete lost particles, if any
OPTION, AUTOPHASE=4;        // Autophase is on, and phase of max energy
                            // gain will be found automatically for cavities
Option, VERSION=10900;

Title, string="Solenoid-1 test";

//----------------------------------------------------------------------------
//Global Parameters

REAL rf_freq             = 1.3e3;     //RF frequency. (Hz)
REAL n_particles         = {n_parts:n};      //Number of particles in simulation.
REAL beam_bunch_charge   = {bunch_charge};      //Charge of bunch. (C)

//Initial Momentum Calculation
REAL Edes    = {energy};  // initial energy in GeV
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

FS_SC: Fieldsolver, TYPE = NONE,
            NX = 16, NY = 16, NZ = 16, // SC grid size is 32^3
            PARFFTX = false,
            PARFFTY = false,
            PARFFTZ = true,  // parallel in the z direction only
            BCFFTX = open,
            BCFFTY = open,
            BCFFTZ = open,
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

TRACK, LINE = SOLine, BEAM = BEAM1, MAXSTEPS = {steps:n}, DT = {dt};
RUN, METHOD = "PARALLEL", BEAM = BEAM1,
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

    return input_file(out_string)
