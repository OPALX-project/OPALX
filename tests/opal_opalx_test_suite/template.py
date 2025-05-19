def get_opal_string(amount:"1e4", steps:"10"):
    return f"""
OPTION, PSDUMPFREQ = 1;   // 6d data written every 300000 time steps (h5).
OPTION, STATDUMPFREQ = 1;  // Beam Stats written every 10 time steps (stat).
OPTION, BOUNDPDESTROYFQ=10; // Delete lost particles, if any
OPTION, AUTOPHASE=4;        // Autophase is on, and phase of max energy
                            // gain will be found automatically for cavities
Option, VERSION=10900;

Title, string="Solenoid-1 test";

//----------------------------------------------------------------------------
//Global Parameters

REAL rf_freq             = 1.3e3;     //RF frequency. (Hz)
REAL n_particles         = {amount};      //Number of particles in simulation.
REAL beam_bunch_charge   = 1e-9;      //Charge of bunch. (C)

//Initial Momentum Calculation
REAL Edes    = 1e-3;  // initial energy in GeV
REAL gamma   = (Edes+EMASS)/EMASS;
REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;    //inital z momentum

//Printing initial energy and momentum to terminal output.
value , {{Edes, P0, OPALVERSION}};

//----------------------------------------------------------------------------
// Solenoids
//
// L:           Physcial element length (m)
// ELEMEDGE:    Physcial start of element (m)
// KS:          Solenoid strength (Rad/m)
// FMAPFM:      Field file (string)

// Note: OPAL scales the field file based on the max magnetic
// field value in the file, not Bz on axis. The max field
// value is normalized to 1 [T], and scaled with KS.
// i.e. The max value in the BF_559 file = 0.162544398 [T].
// Therefore, setting KS = 0.162544398 runs the magnet at max current.

REAL KSBF = 0.162544398;
if (OPALVERSION>10500)
   KSBF = KSBF/1.3528;

BF: Solenoid, L = 0.5, ELEMEDGE=0.0, KS = KSBF, FMAPFN = "BF_550.T7";

SOLine: Line = (BF);

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

TRACK, LINE = SOLine, BEAM = BEAM1, MAXSTEPS = {steps}, DT = {{1.0e-10}};
RUN, METHOD = "PARALLEL-T", BEAM = BEAM1,
    FIELDSOLVER = FS_SC, DISTRIBUTION = Dist;
ENDTRACK;
Quit;

    """

def get_opalx_string(amount="1e4", steps="10"):
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
REAL n_particles         = {amount};      //Number of particles in simulation.
REAL beam_bunch_charge   = 1e-9;      //Charge of bunch. (C)

//Initial Momentum Calculation
REAL Edes    = 1e-3;  // initial energy in GeV
REAL gamma   = (Edes+EMASS)/EMASS;
REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;    //inital z momentum

//Printing initial energy and momentum to terminal output.
value , {{Edes, P0, OPALVERSION}};

//----------------------------------------------------------------------------
// Solenoids
//
// L:           Physcial element length (m)
// ELEMEDGE:    Physcial start of element (m)
// KS:          Solenoid strength (Rad/m)
// FMAPFM:      Field file (string)

// Note: OPAL scales the field file based on the max magnetic
// field value in the file, not Bz on axis. The max field
// value is normalized to 1 [T], and scaled with KS.
// i.e. The max value in the BF_559 file = 0.162544398 [T].
// Therefore, setting KS = 0.162544398 runs the magnet at max current.

REAL KSBF = 0.162544398;
if (OPALVERSION>10500)
   KSBF = KSBF/1.3528;

BF: Solenoid, L = 0.5, ELEMEDGE=0.0, KS = KSBF, FMAPFN = "../BF_550.T7";

SOLine: Line = (BF);

Dist: DISTRIBUTION, TYPE = GAUSS,
        SIGMAX = 0.001,
        SIGMAY = 0.001,
	    SIGMAZ = 0.001;

//----------------------------------------------------------------------------
// Define Field solvers
// The mesh sizes should be a factor of 2
// for most efficient space charge calculation.

FS1: Fieldsolver, Type=NONE,
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

TRACK, LINE = SOLine, BEAM = BEAM1, MAXSTEPS = {steps}, DT = {{1.0e-10}};
RUN, METHOD = "PARALLEL", BEAM = BEAM1,
    FIELDSOLVER = FS1, DISTRIBUTION = Dist;
ENDTRACK;
Quit;
    """


def get_slurm_string(executable, filename):
    return f"""#!/bin/bash
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1               # No. of nodes
#SBATCH --ntasks-per-node=1     # No. of MPI ranks per node. Merlin CPU nodes have 44 cores
#SBATCH --cpus-per-task=1      # No. of OMP threads
#SBATCH --time=00:10:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=nomultithread    # Without hyperthreading
##SBATCH --exclusive            # The allocations will be exclusive if turned on (remove extra hashtag to turn on)

#SBATCH --output={filename}.out  # Name of output file
#SBATCH --error={filename}.err    # Name of error file

export OMP\_NUM\_THREADS=1
export OMP\_PROC\_BIND=spread
export OMP\_PLACES=threads


srun --cpus-per-task=1 {executable} {filename} --info 10
"""
