ReadMe: How to use and edit the SlicedBunch in OPAL


version: 26.08.2008
by: Frank Hantschel

This readme has two main sections: The first provides information how to use OPAL, the second gives more specific information how to extend the existing code. Text in curly brackets {} was extracted from the code itself or the inputfile.  




1) How to use OPAL



For the following section I recommend to use a Linux-shell. To use the sliced bunch, first copy the OPAL directory from home2/hantschel/svnwork/OPAL to your home directory. Then compile OPAL using the "make" command in OPAL/src. Afterwards change to the following directory:

OPAL/opal-Tests/RegressionTests/envelope-test1/OPAL

Then run OPAL typing the following command:

./run

This calls OPAL using the input file "envelope-test1.in". OPAL executes the following steps:

a) initialization of the bunch
b) build up field list and calculate the electric and magnetic fields
c) loop over timesteps
d) write output to "Bop.dat" and "Bop.sli"



a) initialization of the bunch:



OPAL reads in "envelope-test1.in". The following input parameters can be adjusted there: 

- properties of the distribution:

{Dist1:DISTRIBUTION, DISTRIBUTION=gauss,
sigmax=  1.0e-03, sigmapx=1.0e-4, corrx=0.5,
sigmay=  2.0e-03, sigmapy=1.0e-04, corry=0.5, 
sigmat=  3.0e-03, sigmapt=1.0e-4, corrt=0.0, pt=1.23, nbin=10;}

sigmax and sigmay determine the initial width of the bunch in x and y direction, respectively ([m])
sigmat determines the length of the bunch ([m]).
nbin determines the number of slices used (output just for nbin-2 slices!).

All other variables of the distribution are not taken yet.

- properties of the beam:

{beam1: BEAM, PARTICLE=ELECTRON, GAMMA=gamma, pc=P0, NPART=2, BFREQ=1e11, BCURRENT=0.001, CHARGE=-1;}

GAMMA determines the value of gamma
pc determines the initial momentum and therefore the energy as well 
BFREQ,BCURRENT and CHARGE determine the initial charge Q (in C) via 'CHARGE * BCURRENT / BFREQ'. Here, CHARGE has unit of multiples of the unit charge e, BCURRENT has [A] and BFREQ has unit [Hz].

The remaining input parameters are fixed in the source code (default value is given as well):

- TrackRun.cpp (L194-195):

{Bz0=0;}	initial value of the magnetic field 

{cen = 0;}	center of the bunch ([m]).

- Distribution.cpp (L220-222):

{frac=0.0;}	slope parameter (gradient), for rectangle profil 0 

iniBetBunch(..) function:

{mX=0;}		emittance for x coordinate, BET has default value 0

{mY=0;}		emittance for y coordinate, BET has default value 0

SLPartBunch.cpp (L89, 93, 108):

{wfraction=1;}	determines the percentage of output (0: no output, 1: whole bunch is written to file)
 
{i2=0;}		form of the bunch: 0 rectangle, 1 gaussian

{tempB->setSolver(12);} 
tracking method (OnAxis, Radial...), default is 12 (track all). For more information concerning individual modes please consult the source code of the BET programme (system.C and bunch.C)


All these parameter determine the bunch and its properties. To track the bunch through the beamline, OPAL reads in the following to variables from "envelope-test1.in":

{MAXSTEPS=290}	     number of timesteps executed by OPAL
 
{DT=1.0e-11}	     magnitude of one timestep in seconds. Please be aware that the final timestep used can differ from DT, due to dynamic step intergration


Moreover OPAL reads in beamline elements from "envelope-test1.in":

- Solenoid:

{FINEG_MSL10: Solenoid, L=0.3, KS=0.2, FMAPFN="FINEG-MSL10.T7", ELEMEDGE=0.5;}

legend:

{FINEG_MSL10} name of the element (just one element for a certain name)

{Solenoid}    type

{L}	      length (in m), however the used length is already determined by the fieldmap
 
{KS}	      maximum magnetic field strength of Bz in tesla

{FMAPFN}      fieldmap (filename)

{ELEMEDGE}    position (Head!) in [m]


- RFCavity (=SWA!):

{FIRF_01: RFCavity, L=0.018, VOLT=100, FMAPFN="FINLB01-RACF.T7", ELEMEDGE=0.20, TYPE="STANDING", FREQ=1498.956, LAG=0.0;}

legend:

{FIRF_01}     name

{RFCavity}    type

{L}	      magnitude of one timestep in seconds. Please be aware, that the final timestep used can differ from DT, due to dynamic step intergration

{VOLT}	      maximum electric field strength of Ez in [MVolt / m]

{FMAPFN}      fieldmap (filename)

{ELEMEDGE}    position (Head!)

{Type}	      type ("Standing" for SWA)

{FREQ}	      frequency [MHz]

{LAG}         phase in Rad. BEWARE: The phase in BET is always relative to the zero point and therefore includes a shift (w0*z0/c). This is different from OPAL, which defines the phase relative to the starting point of the field.

Until now, only these two types of elements were tested. For the other types the necessary function are already implemented, but these can still contain errors. 



b) build up field list:



OPAL initializes all elements using the input file "envelope-test1.in" and calculates all fields anew every timestep. The fieldstrength itself is written into "Bop.dat" and can be extracted by the analyzation tool SDDS or by gnuplot (Efeld via 'plot "Bop.dat" u 1:30 w l'). More information about the output can be found in section d) of this readme. 



c) loop over timesteps:



As mentioned before, OPAL executes every timestep using the input variable DT. The length of every individual timestep can differ from DF due to dynamic step integration. To track a certain length, one has to assign certain values to DT and MAXSTEPS. To obtain these values, I would recommend to try approximate values for DT and MAXSTEPS (eg using DT*MAXSTEPS=LENGTH) and to determine the exact parameters by try-and-error. After each timestep, OPAL writes output on the screen and to the two files 'Bop.dat' and 'Bop.sli'. OPAL will finish the tracking by the command line " OPAL> execute complete". 



d) output:



The sliced bunch in OPAL creates two output files: Bop.dat and Bop.sli (Bop is an abbreviation for "Bet output"). To analyze the output data, one can use the tool SDDS. To execute SDDS based on OPAL output, use one of the following two command lines. In order to use the SDDS output, add the following command line to your .bashrc file:

export PATH=.:${HOME}/bin:/afs/psi.ch/project/fel/opt/epics/extensions/bin/linux-x86-64:${PATH}

The first plot command is as follows (do not forget to type Bop as input parameter!):

plotAll.sh Bop

If this command does not work, you have to copy the plotAll.sh file into your /bin directory. You will find plotAll.sh in the following directory:

home2/hantschel/bin/plotAll.sh

SDDS creates ten different plots where all of these give values over the whole beamline tracked. One can use the keys "n" (next) and "p" (previous to navigate between individual plots. Please do not use "b", since this command extends the plot over the whole width of the screen. All averages, maxima and minima are determined using the whole bunch length / all different slices.  

- E		       energy of the bunch
- Rx, MaxRx, MinRx     bunch width in x direction with average, minimum and maximum value
- Ry, MaxRy, MinRy     bunch width in y direction with average, minimum and maximum value
- Imax, Irms	       maximum and average value of the current
- tau		       bunch length (unit of time!)
- emtnx, emtny	       emittance for x and y coordinates
- Px, Py	       transversal momentum, hence:  m*dR/dt
- EField	       Efield (z component)
- BField	       BField (z component)
- x0, y0	       OffAxis components (while using OnAxis: x0 = y0 = 0)

The second plot command is as follows:

sddsshow Bop.sli

This gives the profil of the bunch (separates slices). Using the interface, one can select certain parameters and plot them via "Strg-p". I strongly recommmend to enable the options "split pages" and "y-separate" using the plot window, since otherwise different plots interfere with each other. As before, using "n" (next) und "p" (previous) allows to navigate between individual plots.

Additionaly one can use gnuplot to process output data. On the one hand the resulting plots are much more exact and gnuplots allows the user to set x ranges and other plot properties. On the other hand, this is more complicated than to use SDDS: To plot a certain variable, look at the header of the particular output file (Bop.dat or Bop.sli) and get the number of the according column in this file. For example: To get the electric field strength, one would find the variable EField in the 30th column of "Bop.dat". To plot this value, one executes the following command line in gnuplot:   

'plot "Bop.dat" u 1:30 w l'



2) How to extend BET



The following files contain information for the sliced bunch routine:

In svnwork/OPAL/src:

/Track/TrackRun.cpp
/Distribution/Distribution.cpp .h
/Algorithms/ParallelSliceTracker.cpp .h
/Algorithms/bet/bunch.cpp .h
/Algorithms/bet/slice.cpp .h
/Algorithms/bet/SLPartBunch.cpp .h
/Algorithms/bet/physconst.h
/Algorithms/bet/error.cpp .h
/Algorithms/bet/math/*
/Algorithms/bet/libprf/*

In svnwork/OPAL/classic/5.0/src:

/AbsBeamline/Solenoid.cpp .h
/AbsBeamline/Multipole.cpp .h
/AbsBeamline/RBend.cpp .h
/AbsBeamline/RFCavity.cpp .h
/AbsBeamline/TravelingWave.cpp .h
/AbsBeamline/BeamBeam.cpp .h
/AbsBeamline/PartBunch.cpp .h

Some functions of the latter files should not work anymore for the ParticleBunch, to use them one would have to change names and to copy the old function back to files.  

This section is again divided in four subsections:

a) general structure and data
b) initialization
c) timestep 
d) output

a) general structure and data

The original BET programme is based on the files bunch.h and bunch.cpp: They contain the main functions and all data of one sliced bunch. Each bunch is divided into several slices. Hence, the bunch.h includes an array of (pointers of) slices. The following slice parameters are saved in slice.h:

   beta normalized velocity (total) [-]
   z    slice position [m]
   x    beam size x (rms) [m]   
   y    beam size y (rms) [m]
   px   beam divergence x [rad]
   py   beam divergence y [rad]
   X0   position centroid x [m]
   Y0   position centroid y [m]
   pX0  angular deflection centriod x
   pY0  angular deflection centroid y

To access these variables, use s->p[SLI_xxx] in bunch.cpp, where xxx is determined by the following list:

   #define SLI_z     0
   #define SLI_beta  1
   #define SLI_x     2
   #define SLI_px    3
   #define SLI_y     4
   #define SLI_py    5
   #define SLI_x0    6
   #define SLI_px0   7
   #define SLI_y0    8
   #define SLI_py0   9

To access this information outside bunch.cpp, I used public functions of the bunch class in SLPartBunch. 

BEWARE: I did not initialize the R[]() vector as used in the PartBunch class but used get and set functions in SLPartBunch.cpp and bunch.cpp for this purpose.

The bunch class is a member of the class SLPartBunch, which is a subclass of PartBunch. SLPartBunch just calls particular functions of bunch.cpp. Hence, SLPartBunch mainly is used as a shell for the BET bunch class as implemented in bunch.cpp. To call these functions in ParallelSliceTracker.cpp, I also implemented a bunch of virtual functions in the basic class PartBunch.h. The whole process of initialization and executing timesteps is initiated and supervised by TrackRun.cpp. 


b) Initialization


As mentioned before, TrackRun.cpp starts the initialization by calling Distribution::createSlicedBunch() (Distribution.cpp). This function uses the beam parameter as input variables and adds the distribution parameter inside. For a detailed description of these input variables look at the user section of this readme. After gathering all necessary numbers for the initialization, Distribution::createSlicedBunch() calls SLPartBunch::iniBetBunch of SLPartBunch.cpp. This function creates a sliced BET bunch and assignes the input values to this element. To do that I used the following functions of bunch.cpp:     

  setCharge(Q)	     		  [sets charge]

  setEnergy(energy)		  [sets inititialisation energy]

  setEy(0)			  [sets emmitance]

  setLShape(i2?bsGauss:bsRect,center,width,frac) [sets longitudinal shape]

  setTShape(mX,mY,bX,bY,Bz0)  			 [sets transversal shape]

All these functions were not changed in function but copied from the BET programme. All conversions were done before, e.g. to convert the energy from GeV (OPAL input) to eV (as used in BET). For more detailed information concerning these conversions please consult the doxygen documentation. After calling iniBetBunch, the initialization is complete. To obtain changes, I would recommend to edit the function iniBetBunch and not the functions inside bunch.cpp in order to preserve the structure of the programme.


c) timestep


The timestep was mainly implemented using the execute() function of ParallelSliceTracker.cpp. This function builds up a field list and calculates all exterior values for the timestep such as the magnetic and electric field strength. ParallelSliceTracker::execute loops over timesteps and for each calls the function Bunch::run(). The structure of Bunch::run was a again copied from BET but I obtained some changes:

- run() does not build up the field list anymore, this task is done by execute(). 

- the exterior values (Bfield, Efield, K-values) are saved as members of the bunch class (bunch.h), I used an array of vector_t. Bunch::run() gets these values by access on these variables instead of calculating them. 

- BET calculates only 80 % of the bunch, due to special behaviour of slices at the head and tail of the bunch. OPAL does 100 %.  

For more details, take a look at the doxygen documentation of bunch.cpp and ParallelSliceTracker.cpp. Bunch::run() automatically assignes the new values for all slice parameter. 

An important point is the time. In theory, there are three different variables which represent the time: One time in the bunch.cpp, another in ParallelSliceTracker.cpp and the last in TrackRun.cpp. I created functions to update the time. Since the field list is build up differently in OPAL and BET, the fields are not calculated at the same time for the same beamline and the same variables in OPAL and BET. 


d) output


The following three functions did the output in BET:

- bunch::write() -> saves bunch data (single output)

- bunch::writeStats() -> writes into .dat file for each timestep

- bunch::writeSlice() -> writes into .sli file for each timestep

The SDDS-output uses .dat and .sli files, hence, the last two functions are important in order to produce SDDS output. SLPartBunch includes four output functions but uses only one (SLPartBunch::BetOut()). As before, to edit the output I would recommend to edit functions in SLPartBunch.cpp or to create new functions in bunch.cpp. 

For field data (Bfield and Efield) BET uses vectors of the new type "field". Since OPAL uses the type "vector_t", I changed the main variables to the vector_t type. I would recommend to use this type for further implementation.
 
The name of the output files (Bop.dat and Bop.sli) is manually set in TrackRun.cpp (L123). Moreover TrackRun.cpp opens and closes these files, whereas the functions in bunch.cpp take the open file as an input parameter and add data. The header of these files is written at the first call of ParallelSliceTracker::execute() by the corrresponding function above (Bunch::writeStats() for .dat and Bunch::writeSlice() for .sli files). To write field data I created new functions (Bunch::AvBField() and Bunch::AvEField()) to calculate the average B and E field over all slices. For more details consult the doxygen documentation of the three output functions above.   
