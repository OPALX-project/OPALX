/* bunch.h
   bunch class definition

   Includes most of the slice dynamics.
   References to equations relate to:
   HOMDYN STUDY FOR THE LCLS RF PHOTO-INJECTOR
   M. Ferrario, J.E. Clendenin, D.T. Palmer, J.B. Rosenzweig, L. Serafini
   LNF-00/004 and SLAC-PUB 8400 (March 3, 2000)

Project: Beam Envelope Tracker (BET)

Revision history
Date          Description                                     Programmer
------------  --------------------------------------------    --------------
07-03-06      Created                                         Rene Bakker

Last Revision:
$Id: bunch.h 218 2008-03-15 16:01:20Z bakker $


Extended and embedded in OPAL: AMAS

*/


#ifndef _BUNCH_DEF
#define _BUNCH_DEF

#include <stdio.h>
#include "Algorithms/bet/EnvelopeSlice.h"
#include "Algorithms/bet/profile.h"
#include <vector>

enum EnvelopeBunchParameter {
    sp_beta,      /// normalized velocity (total) [-]
    sp_gamma,     /// Lorenz factor
    sp_z,         /// slice position [m]
    sp_I,         /// current [A]
    sp_Rx,        /// beam size x [m]   
    sp_Ry,        /// beam size y [m]
    sp_Px,        /// beam divergence x
    sp_Py,        /// beam divergence y
    sp_x0,        /// position centroid x [m]
    sp_y0,        /// position centroid y [m]
    sp_px0,       /// angular deflection centriod x
    sp_py0        /// angular deflection centroid y
};

// bit flags
// ==================================================================================

enum SolverParameters {
    sv_fixedStep    = 0x0001,  /// solve DV with fixed time-step 
    sv_fieldOutside = 0x0002,  /// solve field outside of DV
    sv_radial       = 0x0004,  /// include radial movement in DV
    sv_offaxis      = 0x0008,  /// include off-axis movement in DV
    sv_lwakes       = 0x0010,  /// longitudinal wakes
    sv_twakes       = 0x0020,  /// transverse wakes
    sv_s_path       = 0x0100   /// track along s-axis (instead of z)
};

enum DataStatus {
    ds_fieldsSynchronized = 0x0001,    /// fields synchronized with other MPI nodes
    ds_slicesSynchronized = 0x0002,    /// slice data synchronized with other MPI nodes
    ds_currentCalculated  = 0x0004,    /// wakes and space-charge fields calculated
    ds_wakesZCalculated   = 0x0008,    /// longitudinal wakes calculated
    ds_wakesXYCalculated  = 0x0010,    /// transverse wakes calculated
    ds_spaceCharge        = 0x001c     /// previos 3 combined (OR)
};

// pre-defined bunch shapes
enum EnvelopeBunchShape {
    bsRect,
    bsGauss
};

// output formats for data
enum OutputFileFormat {
    oFormat_tab,
    oFormat_sdds
};

class EnvelopeBunch {

    /// see enum SolverParameters
    int solver;
    /// see enum DataStatus
    int dStat;
    /// number of slices in bunch
    int n;

    /// local time in bunch (s)
    double t;
    /// accumulated time offset by tReset function
    double t_offset;
    /// intrinsic normalized emittance of slice (m rad)
    double emtnx0,emtny0;
    /// intrinsic normalized emittance Bush effect (m rad)
    double emtbx0,emtby0;
    /// offset of the coordinate system when tracking along the s-axis (m)
    double dx0,dy0;
    /// rotation of coordinate system when tracking along the s-axis (rad)
    double dfi_x,dfi_y;
    /// magnetic field on cathode (T)
    double Bz0;
    /// total bunch charge (C)
    double Q;
    /// average current on creation of bunch (see setLshape)
    double I0avg;


private:
    /* Variables and functions for Integration routine */

    /// allocate memory for an array of n doubles
    double *vector(int);
    /// allocate memory for an array of n integers
    int *ivector(int);


    /// write column header 
    /// file pointer
    /// output format
    /// name of column
    /// type of column (long, double, .....)
    /// units of column
    /// description of column
    void defineColumn(FILE *, OutputFileFormat, char *, char *, char *, char *);              

    /**
      // first element to act on bunch (others auto-linked)
      Element *fE;
    */

    /** current Slice
            - set in run() & cSpaceCharge()
            - used in derivs(); & zcsI();
    */ 
    int cS;               
    
    /// cathode position          
    double zCat;
    /// transverse wake field x
    double *Exw;
    /// transverse wake field y
    double *Eyw;
    /// longitudinal wake field
    double *Ezw;
    /// Longitudinal Space-charge field
    double *Esct;
    /// Transverse Space-charge term:   Eq.(9)
    double *G;


public:
    /// electric field
    Vector_t Esl;
    /// magnetic field
    Vector_t Bsl;
    /// radial focussing term beam
    Vector_t KRsl;
    /// transverse kick of beam 
    Vector_t KTsl;
    
    ///define value of K for each slice (array)
    Vector_t* KR;
    Vector_t* KT; 
    ///define values for B and E field
    Vector_t* EF;
    Vector_t* BF;


    /// current profile of bunch (fit)
    Profile  *I;
    /// array of slices
    EnvelopeSlice **s;
    /// gives the sign of Q
    int sign;


#ifdef USE_MPI
    ///data buffer to synchronize (replaces local storage of slice data) 
    double **mpiBuffer;
    /// size of mpi buffers
    int mpiBsize;

    /// relocate slice parameter memory
    void relocate();
    /// synch slices with other mpi processes
    void sync_slices(int);
    /// sync space-charge & wake fields for output
    void sync_cFields();
#endif

    // run statistics on slices
    // select parameter to run stats on
    // fractional window
    void runStats(EnvelopeBunchParameter sp, double w, double *xAvg, double *xMax, double *xMin, double *sigma, int *nValid);

    // calculate the emittance
    // fractional window
    // emt x (normalized)
    // emt y (normalized)
    // number of valid points
    void calcEmittance(double, double *, double *, int*);

    // calculate the energy chirp and uncorrelated energy spread
    // fractional window
    // average gamma
    // chirp: dgdt
    // incoherent energy spread
    // number of valid points
    void calcEnergyChirp(double, double*, double*, double*, int*);

    void cSpaceCharge(); // Calculate space-charge fields

    void printSC(const char *fName); //print SC data to file

    void writeSC(int slices); // Write SC factors

    void DefineSlices(int no, double beta, double pos, double bsX, double bsY);

    void writeSlices();

    double getGamma(int i);

    double getBeta(int i);

    void setZ(int i, double coo); // set Z coordinate

    double getZ(int i); // get ith Z coordinate

    double getX(int i); // get ith X coordinate

    double getY(int i); // get ith Y coordinate

    double getX0(int i); // get ith X coordinate

    double getY0(int i); // get ith Y coordinate

    double getPx(int i); // get momentum 

    double getPy(int i);

    double getPz(int i);

    double getPx0(int i);

    double getPy0(int i);

    double getT();

    void setEx(double emi); // set value for the intrinsic emittance

    void setEy(double emi);

    void plotR();

    double AvBField();

    double AvEField();


    /// allocate values for slices
    /// calculate the current current distribution
    void calcI();
    /**
      void cWake();        // Calculate wake fields 
    */
public:

    /// create a bunch with n slices
    EnvelopeBunch(int = 101);
    /// read a bunch distribution from file
    EnvelopeBunch(char *);

    // read a bunch distribution from a homdyn HSAVE.COR file
    // filename
    // bunch charge (C)
    // slice emittance x (m rad)
    // slice emittance y (m rad)
    // magnetic field on cathode (T)
    EnvelopeBunch(char *, double, double, double, double);
    ~EnvelopeBunch();

    // information 
    /// the number of slices in the bunch
    int getN();


    // diagnostic functions (1)
    /// calculate the average energy of the bunch
    double Eavg();

    // diagnostic functions (2)
    /// calculate <z> [m]
    double zAvg();
    /// calculate tail of bunch [m]
    double zTail();
    /// calculate the head of the bunch [m]
    double zHead();
    /// read time-stamp of bunch
    double time();
    
    /// read bunch from file
    void read(char *);

    /// write bunch to file (NULL - stdout)
    void write(const char * = NULL);
    /// write to filepointer
    void write(FILE *);

    // write EnvelopeBunch statistics 
    // file pointer
    // fractional window for statistics
    void writeStats(FILE *, double, OutputFileFormat);
    
    // write all slice info to file
    // file pointer
    // format of output file
    void writeSlice(FILE *, OutputFileFormat); 

    // property defining functions

    /// set the charge of the bunch
    void setCharge(double);

    // set longitudinal shape of bunch
    // shape of bunch (rect/gauss)
    // z [m] of center 
    /* pos. length [m], neg. duration [s]
       (rms value) */
    // for Gauss only: fraction of Gauss used
    void setLShape(EnvelopeBunchShape, double, double, double = 0.95);                    

    // set transverse shape of bunch
    // normalized emittance x [m rad]
    // normalized emittance y [m rad]
    // Rx  [m]
    // Ry  [m]
    // Bz0 [T]
    void setTShape(double, double, double, double, double);          

    // set transverse offset of bunch
    // x0  [m]
    // px0 [rad]
    // y0  [m]
    // py0 [rad]
    void setTOffset(double, double, double, double);                 

    // set particle energy of bunch
    // energy [eV]
    // correlated energy spread [eV/m]
    void setEnergy(double, double = 0.0);                     
    
    /// set the DE solver flag
    void setSolver(int);

    /// move the complete bunch forward such that the head of the bunch matches the cahtode position
    double moveZ0(double);

    /// backup slice values
    void backup();

    // check if present position is marked as a screen, returns 1 if written
    // present position
    // file pointer for output
    // format of output file
    int checkScreen(double, FILE *, OutputFileFormat);                

    /// reset time of bunch (returns the offset applied)
    /// time difference (0.0 - auto-sync)
    double tReset(double = 0.0);

    // solving of the differential equations
    // RH equation 
    /// t, Y[], dYdt[]
    void derivs(double, double *, double *);

    /// update space-charge fields
    void updateFields();

    /// do the calculation 
    /// zCat
    void run(double, double);  

    //void runSS(double, double);

};

typedef EnvelopeBunch *EnvelopeBunchP;

#endif
