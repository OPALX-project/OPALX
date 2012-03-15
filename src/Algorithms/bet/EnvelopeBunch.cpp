/* bunch.C
   bunch class definition

Project: Beam Envelope Tracker (BET)

Revision history
Date          Description                                     Programmer
------------  --------------------------------------------    --------------
07-03-06      Created     




Rene Bakker

Last Revision:
$Id: bunch.C 285 2008-05-30 13:23:26Z bakker $
*/

#define SVN_DATE "$Date: 2008-05-30 15:23:26 +0200 (Fri, 30 May 2008) $"

#ifdef USE_MPI
#include"mympi.h"
#endif

#include "Ippl.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <Physics/Physics.h>

#include <vector>
#include "Solvers/FFTPoissonSolver.h"

#include "Algorithms/bet/error.h"             // message handling
#include "Algorithms/bet/math/functions.h"    // special functions
#include "Algorithms/bet/math/root.h"         // root finding routines
#include "Algorithms/bet/math/sort.h"         // sorting routines
#include "Algorithms/bet/math/integrate.h"    // integration routines
#include "Algorithms/bet/math/linfit.h"       // linear fitting routines
#include "Algorithms/bet/math/savgol.h"       // savgol smoothing routine
#include "Algorithms/bet/math/rk.h"           // Runge-Kutta Integration
#include "Algorithms/bet/EnvelopeBunch.h"             // bunch class definition 

// #include "utils.h"        // general utilities

extern Inform* gmsg;

/* for space charge
   ignore for too low energies
   beta = 0.05 -> 640 eV
   beta = 0.10 ->   3 keV
   beta = 0.20 ->  11 keV
   beta = 0.30 ->  25 keV
   beta = 0.40 ->  47 keV
   beta = 0.50 ->  90 keV
   beta = 0.60 -> 128 keV 
   beta = 0.70 -> 205 keV
   beta = 0.75 -> 261 keV
   beta = 0.80 -> 341 keV
   beta = 0.85 -> 460 keV */

#define BETA_MIN1 0.30     // minimum beta-value for space-charge calculations: start
#define BETA_MIN2 0.45     // minimum beta-value for space-charge calculations: full impact

/* local functions: RK integration
   ======================================================================== */

/** Dirty programming
 *  Allows odeint (rk.C) to be called with a class member function 
 */

static EnvelopeBunch *itsBunch   = NULL;  // pointer to access calling bunch

/// derivs for RK routine
static void Gderivs(double t,double Y[],double dYdt[]) {
    itsBunch->derivs(t,Y,dYdt);
}

static double Gcur(double z) {
    return itsBunch->I->get(z,itype_lin);
}

/* local functions: root finding (single value for Gauss distribution)
   ======================================================================== */

static double rootValue = 0.0;

static void erfRoot(double x,double *fn, double *df)
{
    double v = erfc(fabs(x));
    double eps = 1.0e-05;

    *fn = v - rootValue;
    *df = (erfc(fabs(x)+eps)-v)/eps;
}


/* private functions
   ======================================================================== */
double *EnvelopeBunch::vector(int nn)
{
    double *r;
    r = (double *) malloc(sizeof(double)*nn);
    if (!r) 
        writeError();
    return r;
}

int *EnvelopeBunch::ivector(int nn)
{
    int *r;
    r = (int *) malloc(sizeof(int)*nn);
    if (!r) 
        //   writeError(errModeAll,errFatal,
        //	       "Memory allocation error EnvelopeBunch::ivector(%d)",nn);
        writeError();
    return r;
}

void EnvelopeBunch::defineColumn(FILE *f,OutputFileFormat oFormat, char *nameStr, char *typeStr, char *unitStr, char *descStr) 
{
    switch (oFormat) {
    case oFormat_tab :
        fprintf(f,"%s%s%s \t", nameStr,unitStr?"-":"",unitStr?unitStr:"");
        break;
    case oFormat_sdds :
        fprintf(f,"&column name=%s, type=%s, units=\"%s\", description=\"%s\" &end\n",
                nameStr,typeStr,unitStr?unitStr:"-",descStr?descStr:"-");
        break;
    }
}

/// run statistics on slices
/// select parameter to run stats on
void EnvelopeBunch::runStats(EnvelopeBunchParameter sp, double w, double *xAvg, double *xMax, double *xMin, double *sigma, int *nValid)
{
    int i,i0,i1,nV = 0;
    double M1,M2,min,max,*v = vector(n);

    switch (sp) {
    case sp_beta:      // normalized velocity (total) [-]
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_beta];
        }
        break;
    case sp_gamma:     // Lorenz factor
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->gamma();
        }
        break;
    case sp_z:         // slice position [m]
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_z];
        }
        break;
    case sp_I:         // slice position [m]
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_beta] > BETA_MIN1) && ((s[i]->p[SLI_z] > zCat) && s[i]->valid)) 
                v[nV++] = I->get(s[i]->p[SLI_z],itype_lin);
        }
        break;
    case sp_Rx:        // beam size x [m]   
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = 2.0*s[i]->p[SLI_x];
        }
        break;
    case sp_Ry:        // beam size y [m]
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = 2.0*s[i]->p[SLI_y];
        }
        break;
    case sp_Px:        // beam divergence x
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_px];
        }
        break;
    case sp_Py:        // beam divergence y
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_py];
        }
        break;
    case sp_x0:        // position centroid x [m]
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_x0];
        }
        break;
    case sp_y0:        // position centroid y [m]
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_y0];
        }
        break;
    case sp_px0:       // angular deflection centriod x
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_px0];
        }
        break;
    case sp_py0:      // angular deflection centroid y
        for (i=1; i<n-1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_py0];
        }
        break;
    default :
        //    writeError(errModeAll,errFatal,"EnvelopeBunch::runStats() Undefined label (programming error)");
        writeError();
        break;
    }

    if (nV <= 0) {
        *xAvg   = 0.0; *xMax = 0.0; *xMin = 0.0; *sigma = 0.0;
        *nValid = 0;
    } else {
        i0 = (int) (0.5*(1.0 - w)*nV);
        if (i0 < 0)   i0 = 0;
        if (i0 >= nV) i0 = nV - 1;
        i1 = nV - i0;
        if (i1 < i0) i1  = i0 + 1;
        nV = i1 - i0;

        M1 = v[i0];
        M2 = v[i0]*v[i0];
        min = v[i0]; max = v[i0];
        for (i=i0+1; i<i1; i++) {
            M1 += v[i];
            M2 += v[i]*v[i];
            if (v[i] < min) min = v[i];
            if (v[i] > max) max = v[i];
        }
        *xAvg   = M1/nV;
        *xMax   = max;    *xMin = min;
        *sigma  = sqrt((M2/nV) - pow(M1/nV,2));
        *nValid = nV;
    }
    free(v);
}

/// calculate the emittance
void EnvelopeBunch::calcEmittance(double w, double *emtnx, double *emtny, int *nValid)
{
    int i,i0,i1,nV;
    double pbc,bg,sx,sxp,sxxp,sy,syp,syyp;
    EnvelopeSlice *sp;

    /** find the amount of active slices
    */
    nV = 0;
    for (i=0; i<n; i++) {
        if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) ++nV;
    }

    if (nV > 0) {
        i0 = (int) (0.5*(1.0 - w)*nV);
        if (i0 < 0)   i0 = 0;
        if (i0 >= nV) i0 = nV - 1;
        i1 = nV - i0;
        if (i1 < i0) i1  = i0 + 1;

        nV = 0;
        bg = 0.0;
        sx = 0.0; sxp = 0.0; sxxp = 0.0;
        sy = 0.0; syp = 0.0; syyp = 0.0;
        for (i=i0; i<i1; i++) {
            if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) {
                ++nV;

                if (solver & sv_radial) {
                    sp = s[i];
                    bg = sp->p[SLI_beta]*sp->gamma();

                    pbc   = bg*sp->p[SLI_px]/(sp->p[SLI_beta]*Physics::c);
                    sx   += pow(sp->p[SLI_x],2); 
                    sxp  += pow(pbc,2); 
                    sxxp += (sp->p[SLI_x]*pbc);

                    pbc   = bg*sp->p[SLI_py]/(sp->p[SLI_beta]*Physics::c);
                    sy   += pow(sp->p[SLI_y],2); 
                    syp  += pow(pbc,2); 
                    syyp += (sp->p[SLI_y]*pbc);
                }
            }
        }
    }

    if (nV == 0) {
        *emtnx   = 0.0; *emtny = 0.0; 
        *nValid = 0;
    } else {
        sx = sx/nV; sxp = sxp/nV; sxxp = sxxp/nV;
        sy = sy/nV; syp = syp/nV; syyp = syyp/nV;

        *emtnx = sqrt((sx*sxp) - pow(sxxp,2) + pow(emtnx0,2) + pow(emtbx0,2));
        *emtny = sqrt((sy*syp) - pow(syyp,2) + pow(emtny0,2) + pow(emtby0,2));
        *nValid = nV;
    }
}


/** calculate the energy chirp and uncorrelated energy spread
*/

// fractional window
// average gamma
// chirp: dgdt
// incoherent energy spread
// number of valid points
void EnvelopeBunch::calcEnergyChirp(double w, double *g0, double *dgdt, double *gInc, int *nValid)   
{
    int i,nV;
    double *dt = vector(n);  /// z --> time
    double *b  = vector(n);  /// beta
    double *g  = vector(n);  /// gamma

    // defaults
    *g0   = 1.0;
    *dgdt = 0.0;
    *gInc = 0.0;

    double zAvg = 0.0;
    double gAvg = 0.0;
    int j = 0;
    for (i=0; i<n; i++) {
        EnvelopeSlice *cs = s[i];
        if ((cs->valid) && (cs->p[SLI_z] > zCat)) {
            zAvg    += cs->p[SLI_z];
            dt[i-j]  = cs->p[SLI_z];
            b[i-j]   = cs->p[SLI_beta];
            g[i-j]   = cs->gamma();
            gAvg    += g[i-j];
        } else ++j;
    }
    nV = n - j;

    if (nV > 0) {
        gAvg = gAvg/nV;
        zAvg = zAvg/nV;
        *g0  = gAvg;
    }

    if (nV > 2) {
        int
            i0,i1;
        double
            dum1,dum2,dum3,dum4,
            rms,gZero,gt;

        i0 = (int) (0.5*(1.0 - w)*nV);
        if (i0 < 0)    i0  = 0;
        if (i0 > nV-2) i0 = nV - 2;
        i1 = nV - i0;

        if ((i1-i0) > 2) {
            // convert z to t
            for (i=i0; i<i1; i++) {
                dt[i] = (dt[i] - zAvg)/(b[i]*Physics::c);
            }

            /// chrip and uncorrelated energy sread
            linfit(&dt[i0],&g[i0],i1-i0,&gZero,&gt,&dum2,&dum3,&dum4);
            *dgdt = gt;

            rms = 0.0;
            for (i=i0; i<i1; i++) {
                rms += pow(g[i]-gZero-gt*dt[i],2);
            }
            *gInc = sqrt(rms/(i1-i0));
        }
    }
    free(g); 
    free(b); 
    free(dt);
}


#ifdef USE_MPI
void EnvelopeBunch::relocate() 
{
    double *buf;
    int i,j,k;

    mpiBsize = SLNPAR*((n/mpi_size) + 1);

    mpiBuffer = (double **) malloc(sizeof(double *)*mpi_size);
    if (mpiBuffer == NULL) {
        writeError(errModeAll,errFatal,"EnvelopeBunch::relocate() - Insufficient memory");
    }

    for (i=0; i<mpi_size; i++) {
        /** last elements used to send extra info through bBcast
        */

        buf = vector(mpiBsize+1);

        k = 0;
        for (j=i; j<n; j += mpi_size) {
            s[j]->setBuffer(&buf[k]);
            k += SLNPAR;
        }
        mpiBuffer[i] = buf;
    }
}

void EnvelopeBunch::sync_slices(int cnt) 
{
    int i;
    double dCnt;

    if ((dStat & ds_slicesSynchronized) == 0) {

        //    writeError(errModeAll,errMessage,"EnvelopeBunch::sync_slices(%d) called",cnt);

        dCnt  = (double) cnt;
        mpiBuffer[mpi_rank][mpiBsize] = dCnt;

        /// make sure we are in sync
        /// MPI::COMM_WORLD.Barrier();
        for (i=0; i < mpi_size; i++) {
            try {
                MPI::COMM_WORLD.Bcast(mpiBuffer[i],mpiBsize+1,MPI_DOUBLE,i);
            } catch (int e) {
                writeError(errModeAll,errFatal,
                        "EnvelopeBunch::mpi_sync(%d) MPI Bcast (i = %d) failed (%d)",cnt,i,e);
            }
            if (mpiBuffer[i][mpiBsize] != dCnt) {
                writeError(errModeAll,errFatal,
                        "EnvelopeBunch::mpi_sync(%d) Synchronization lost (i = %d) %e %e",
                        cnt,i,dCnt,mpiBuffer[i][mpiBsize]);
            }
        }
    }
    dStat |= ds_slicesSynchronized;
    dStat &= (!ds_spaceCharge);
    calcI();
}

/** collect all space charge calculations at master
*/

void EnvelopeBunch::sync_cFields() 
{
    double *buf;

    if ((dStat & ds_fieldsSynchronized) == 0) {
        buf = vector(n);

        //    writeError(errModeAll,errMessage,"EnvelopeBunch::sync_cFields() called");

        try {
            MPI::COMM_WORLD.Allreduce(Esct,buf,n,MPI_DOUBLE,MPI_SUM);
            memcpy(Esct,buf,sizeof(double)*n);

            MPI::COMM_WORLD.Allreduce(G,buf,n,MPI_DOUBLE,MPI_SUM);
            memcpy(G,buf,sizeof(double)*n);

            MPI::COMM_WORLD.Allreduce(Ezw,buf,n,MPI_DOUBLE,MPI_SUM);
            memcpy(Ezw,buf,sizeof(double)*n);

            MPI::COMM_WORLD.Allreduce(Exw,buf,n,MPI_DOUBLE,MPI_SUM);
            memcpy(Exw,buf,sizeof(double)*n);

            MPI::COMM_WORLD.Allreduce(Eyw,buf,n,MPI_DOUBLE,MPI_SUM);
            memcpy(Eyw,buf,sizeof(double)*n);
        } catch (int e) {
            writeError(errModeAll,errFatal,
                    "EnvelopeBunch::mpi_sync_cFields() MPI Bcast send failed: %d",e);
        }

        free(buf);
    }
    dStat |= ds_fieldsSynchronized;
}

#endif

/* public functions
   ======================================================================== */

EnvelopeBunch::EnvelopeBunch(int nSlice) 
{
    n = nSlice;

    KR = new Vector_t[n]; //initialise KR and KT
    KT = new Vector_t[n];

    EF = new Vector_t[n]; //initialise EF and BF
    BF = new Vector_t[n];

    //for (int counter=0;counter<n;counter++) {
    for(int i=0; i<n; i++) {
        KR[i] = (0,0,0);
        KT[i] = (0,0,0);
        BF[i] = (0,0,0);
        EF[i] = (0,0,0);
    }

    //  fE     = NULL;                  // no beam elements defined yet

    // set default DE solver method
    solver = sv_radial | sv_offaxis | sv_lwakes | sv_twakes;

    if (n < 10) {
        //    writeError(errModeMaster,errWarning,"EnvelopeBunch::EnvelopeBunch(%d) called with insufficient slices. Set to 101",nSlice);
        writeError();
        n = 101;
    }

    // defaults: 
    Q        = 0.0;                   // no charge
    t        = 0.0;                   // t = 0 s
    t_offset = 0.0;                   // offset time by tReset function
    emtnx0   = 0.0; emtny0 = 0.0;     // no intrinsic emittance
    emtbx0   = 0.0; emtby0 = 0.0;     // no intrinsic emittance Bush effect
    Bz0      = 0.0;                   // no magnetic field on creation (cathode)

    dx0      = 0.0; dy0    = 0.0;     // no offset of coordinate system
    dfi_x    = 0.0; dfi_y  = 0.0;     // no rotation of coordinate system

    s = new EnvelopeSliceP[n];
    for(int i=0; i<n; i++) s[i] = new EnvelopeSlice();

#ifdef USE_MPI
    mpiBsize = n;
    if (mpi_size > 1) relocate();
#endif

    /// allocate memory for internal arrays
    Esct = vector(n);
    G    = vector(n);
    Exw  = vector(n);
    Eyw  = vector(n);
    Ezw  = vector(n);
    for(int i=0; i<n; i++) {
        Esct[i] = 0.0; G[i]   = 0.0;
        Ezw[i]  = 0.0; 
        Exw[i]  = 0.0; Eyw[i] = 0.0;
    }

    I     = NULL;
    I0avg = 0.0;
    dStat = ds_fieldsSynchronized | ds_slicesSynchronized;
}

EnvelopeBunch::EnvelopeBunch(char *fname) 
{
    n = 0;                         // no slices (yet)
    //  fE     = NULL;                      // no beam elements defined yet

    /// set default DE solver method
    solver = sv_radial | sv_offaxis | sv_lwakes | sv_twakes;

    // writeError(errModeMaster,errMessage,"Loading bunch: %s",fname);
    writeError();
    read(fname);

#ifdef USE_MPI
    mpiBsize = n;
    if (mpi_size > 1) relocate();
#endif

    I     = NULL;
    I0avg = 0.0;
    dStat = ds_fieldsSynchronized | ds_slicesSynchronized;
}


/** create bunch from HOMDYN output file
*/
EnvelopeBunch::EnvelopeBunch(char *fname,double _Q,double ex,double ey,double B0) 
{
    FILE *f;
    int i,j,nLine;
    double y[11],z0;
    EnvelopeSlice *cs;

    /// open the input file and determine the number of slices
    f = fopen(fname,"r");
    if (f == NULL) {
        //  writeError(errModeAll,errFatal,"EnvelopeBunch::read(%s) - cannot open input put file",fname);
        writeError();
    }

    if (fscanf(f,"%d %lf",&i,&z0) < 2) {
        //  writeError(errModeAll,errFatal,
        //       "EnvelopeBunch::read(%s) - Not HSAVE.COR format",fname);
        writeError();
    }
    if (i != 1) {
        //writeError(errModeAll,errFatal,
        //       "EnvelopeBunch::read(%s) - HSAVE.COR only one bunch allowed (%d)",
        //       fname,i);
        writeError();
    }
    nLine = 0;
    while (!feof(f) && (fscanf(f,"%lf",&y[0]) == 1)) ++nLine;
    fclose(f);

    /// determine number of slices
    n = (nLine - 6)/10;

    if (n < 10) {
        //writeError(errModeAll,errFatal,
        //       "EnvelopeBunch::EnvelopeBunch(%s) HSAVE.COR file contains insufficient slices: %d < 10",
        //       n);
        writeError();
    }

    // defaults: 
    //  fE     = NULL;                  // no beam elements defined yet

    Q        = _Q; 
    t        = z0/Physics::c;
    t_offset = 0.0;                   // offset time by tReset function

    emtnx0 = ex;  emtny0 = ey;
    emtbx0 = 0.0; emtby0 = 0.0;     // no intrinsic emittance Bush effect
    Bz0    = B0; 

    dx0    = 0.0; dy0    = 0.0;     // no offset of coordinate system
    dfi_x  = 0.0; dfi_y  = 0.0;     // no rotation of coordinate system

    /// set default DE solver method
    solver = sv_radial | sv_offaxis | sv_lwakes | sv_twakes;

    /// define slices
    s = new EnvelopeSliceP[n];
    f = fopen(fname,"r");
    fscanf(f,"%d %lf",&i,&z0);
    for (i=0; i<6; i++) fscanf(f,"%lf",&y[0]);
    for (i=0; i<n; i++) {
        cs = new EnvelopeSlice();
        for (j=0; j<10; j++) fscanf(f,"%lf",&y[j+1]);

        cs->p[SLI_z]     = y[4];
        cs->p[SLI_beta]  = y[3];
        cs->p[SLI_x]     = y[5]/2.0;
        cs->p[SLI_px]    = y[6]/2.0;
        cs->p[SLI_y]     = y[1]/2.0;
        cs->p[SLI_py]    = y[2]/2.0;
        cs->p[SLI_x0]    = y[7];
        cs->p[SLI_px0]   = y[8];
        cs->p[SLI_y0]    = y[9];
        cs->p[SLI_py0]   = y[10];

        s[i] = cs;
    }
    fclose(f);

#ifdef USE_MPI
    mpiBsize = n;
    if (mpi_size > 1) relocate();
#endif

    /// allocate memory for internal arrays
    Esct = vector(n);
    G    = vector(n);
    Exw  = vector(n);
    Eyw  = vector(n);
    Ezw  = vector(n);
    for (i=0; i<n; i++) {
        Esct[i] = 0.0; G[i]   = 0.0;
        Ezw[i]  = 0.0; 
        Exw[i]  = 0.0; Eyw[i] = 0.0;
    }

    I    = NULL;
    I0avg= 0.0;
    dStat = ds_fieldsSynchronized | ds_slicesSynchronized;
}

EnvelopeBunch::~EnvelopeBunch() 
{
    int i;

    for (i=0; i<n; i++) 
        delete s[i];
    delete s;
    n = 0;

    if (Exw)  free(Exw);
    if (Eyw)  free(Eyw);
    if (Ezw)  free(Ezw);
    if (Esct) free(Esct); 
    if (G)    free(G);
    if (I)    delete I;

#ifdef USE_MPI
    for (i=0; i<mpi_size; i++) {
        free(mpiBuffer[i]);
    }
    free(mpiBuffer);
#endif
}

int EnvelopeBunch::getN() 
{
    return n;
}


double EnvelopeBunch::Eavg() 
{
    int i,nValid;
    double sum;

    nValid = 0;
    sum    = 0.0;
    for (i=0; i<n; i++) {
        if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) {
            sum += s[i]->gamma();
            ++nValid;
        }
    }
    return (nValid > 0?((Physics::EMASS*Physics::c*Physics::c/Physics::q_e)*((sum/nValid)-1.0)):0.0);
}

double EnvelopeBunch::zAvg() 
{
    int i,nV;
    double sum;

    nV  = 0;
    sum = 0.0;
    for (i=0; i<n; i++) {
        if (s[i]->valid) {
            sum += s[i]->p[SLI_z];
            ++nV;
        }
    }
    if (nV < 1) {
        //  writeError(errModeAll,errFatal,"EnvelopeBunch::zAvg() no valid slices left n=%d",
        //       n);
        writeError();
    }
    return (sum/nV);
}

double EnvelopeBunch::zTail() 
{
    int i;
    double min;

    i = 0;
    while ((i < n) && (!s[i]->valid)) ++i;
    if (i == n) {
        // writeError(errModeAll,errFatal,"EnvelopeBunch::zTail() no valid slices left n=%d",
        //       n);
        writeError();
    } else {
        min = s[i]->p[SLI_z];
    }
    for (i=i+1; i<n; i++) 
        if ((s[i]->p[SLI_z] < min) && (s[i]->valid)) 
            min = s[i]->p[SLI_z];
    return min;
}

double EnvelopeBunch::zHead() 
{
    int i;
    double max;

    i = 0;
    while ((i < n) && (s[i]->valid == 0)) ++i;
    if (i == n) {
        // writeError(errModeAll,errFatal,"EnvelopeBunch::zTail() no valid slices left");
        writeError();
    } else {
        max = s[i]->p[SLI_z];
    }
    for (i=1; i<n; i++) 
        if (s[i]->p[SLI_z] > max) max = s[i]->p[SLI_z];
    return max;
}

double EnvelopeBunch::time() 
{
    return t;
}

void EnvelopeBunch::read(char *fname) 
{
    FILE *f;
    char line[2048];
    int  i,nArg,lineNr;

    /// delete old slice definition
    if (n > 0) {
        for (i=0; i<n; i++) delete s[i];
        delete s;

        if (Exw)  {free(Exw);  Exw  = NULL;}
        if (Eyw)  {free(Eyw);  Eyw  = NULL;}
        if (Ezw)  {free(Ezw);  Ezw  = NULL;}
        if (Esct) {free(Esct); Esct = NULL;}
        if (G)    {free(G);    G    = NULL;}
        if (I)    {delete I;   I    = NULL;}
    }

    f = fopen(fname,"r");
    lineNr = 0;

    if (f == NULL) {
        // writeError(errModeAll,errFatal,"EnvelopeBunch::read(%s) - cannot open input put file",fname);
        writeError();
    }

    fgets(line,2048,f); ++lineNr;
    switch (sscanf(line,"%lf %lf %*s\n",&t,&t_offset)) {
    default :
        break;
    case 1 :
        t_offset = 0.0; // old bunch file has no offset info
        break;
    case 0:
        //writeError(errModeAll,errFatal,
        // "EnvelopeBunch::read(%s) - Cannot read time-stamp on line nr %d",
        //  fname,lineNr);
        writeError();
        break;
    }

    fgets(line,2048,f); ++lineNr;
    nArg = sscanf(line,"%d %*s\n",&n);
    if ((nArg < 1) || (n < 1)) {
        // writeError(errModeAll,errFatal,
        //"EnvelopeBunch::read(%s) - Cannot read numer of slices on line nr %d",
        // fname,lineNr);
        writeError();
    }

    fgets(line,2048,f); ++lineNr;
    nArg = sscanf(line,"%lf %*s\n",&Q);
    if ((nArg < 1) || (Q < 0.0)) {
        //writeError(errModeAll,errFatal,
        //"EnvelopeBunch::read(%s) - Illegal EnvelopeBunch Charge (Q=%lg C) on line nr %d",
        // fname,Q,lineNr);
        writeError();
    }

    fgets(line,2048,f); ++lineNr;
    nArg = sscanf(line,"%lf %lf %*s\n",&emtnx0,&emtny0);
    if ((nArg < 2) || (emtnx0 < 0.0) || (emtny0 < 0.0)) {
        // writeError(errModeAll,errFatal,
        //	       "EnvelopeBunch::read(%s) - Illegal intrinsic emittance defined: (emtn0=[%lg,%lg] m rad) on line nr %d",
        //     fname,emtnx0,emtny0,lineNr);
        writeError();
    }

    fgets(line,2048,f); ++lineNr;
    nArg = sscanf(line,"%lf %lf %*s\n",&emtbx0,&emtby0);
    if ((nArg < 2) || (emtbx0 < 0.0) || (emtby0 < 0.0)) {
        // writeError(errModeAll,errFatal,
        //	       "EnvelopeBunch::read(%s) - Illegal intrinsic Bush emittance defined: (emtn0=[%lg,%lg] m rad) on line nr %d",
        //	       fname,emtbx0,emtby0,lineNr);
        writeError();
    }

    fgets(line,2048,f); ++lineNr;
    nArg = sscanf(line,"%lf %*s\n",&Bz0);
    if (nArg < 1) {
        //writeError(errModeAll,errFatal,
        //	       "EnvelopeBunch::read(%s) - Illegal Magnetic Field (Bz0=%lg T) on line nr %d",
        //	       fname,Bz0,lineNr);
        writeError();
    }

    fgets(line,2048,f); ++lineNr;
    nArg = sscanf(line,"%lf %lf %*s\n",&dx0,&dy0);
    if (nArg < 2) {
        //writeError(errModeAll,errFatal,
        //	       "EnvelopeBunch::read(%s) - Error reading average bunch offset: (x0=%lg m, y=%lg m) on line nr %d",
        //	       fname,dx0,dy0,lineNr);
        writeError();
    }

    fgets(line,2048,f); ++lineNr;
    nArg = sscanf(line,"%lf %lf %*s\n",&dfi_x,&dfi_y);
    if (nArg < 2) {
        // writeError(errModeAll,errFatal,
        //	       "EnvelopeBunch::read(%s) - Error reading average bunch rotation: (x0=%lg m, y=%lg m) on line nr %d",
        //   fname,dfi_x,dfi_y,lineNr);
        writeError();
    }

    fgets(line,2048,f); ++lineNr; // header of slices

    /// read all slices
    s = new EnvelopeSliceP[n];
    for (i=0; i<n; i++) {
        int nScan;

        ++lineNr;
        s[i]  = new EnvelopeSlice();
        nScan = s[i]->read(f);

        if (nScan < (SLNPAR+1)) {
            // writeError(errModeAll,errFatal,"EnvelopeBunch::read(%s) - Cannot read slice %d on line nr %d (%d args)",fname,i,lineNr,nScan);
            writeError();  
        }
    }
    fclose(f);

    Exw = vector(n);
    Eyw = vector(n);
    Ezw = vector(n);
    Esct= vector(n);
    G   = vector(n);
    for (i=0; i<n; i++) {
        Esct[i] = 0.0; G[i]   = 0.0;
        Ezw[i]  =  0.0; 
        Exw[i]  = 0.0;  Eyw[i] = 0.0;
    }

    I = NULL;
}

void EnvelopeBunch::write(const char *fname) 
{
    FILE *f;

    f = fname?fopen(fname,"w"):stdout;

    if (!f) {
        // writeError(errModeAll,errGeneral,"EnvelopeBunch::write(%s) - cannot open output file",fname);
        writeError();
        f = stdout;
    }
    write(f);

    if (f != stdout) fclose(f);
}

void EnvelopeBunch::write(FILE *f) 
{
    int i,j,nV;

    fprintf(f,"%20.14le \t %20.14le \t t \t Time stamp and offset time of bunch\n",
            t,t_offset);
    fprintf(f,"%20d \t N \t Number of slices\n",n);
    fprintf(f,"%20.14le \t Q \t EnvelopeBunch Charge [C]\n",Q);
    fprintf(f,"%20.14le \t %20.14le \t emtn0 \t Intrinsic slice emittance (x,y) [m rad]\n",
            emtnx0,emtny0);
    fprintf(f,"%20.14le \t %20.14le \t emtn0 \t Intrinsic slice emittance Bush effect (x,y) [m rad]\n",
            emtbx0,emtby0);
    fprintf(f,"%20.14le \t Bz0 \t Magnetic field on creation (cathode) [T]\n",Bz0);
    fprintf(f,"%20.14le \t %20.14le \t x0,y0 \t coordinate offset   in s-axis mode [m]\n",
            dx0,dy0);
    fprintf(f,"%20.14le \t %20.14le \t fi_x,fi_y \t coordinate rotation in s-axis mode [deg]\n",
            180.0*dfi_x/Physics::pi,180.0*dfi_y/Physics::pi);
    fprintf(f,"#"); s[0]->writeHeader(f);
    for (i=0; i<n; i++) 
        s[i]->write(f);
}

//FIXME: /** ---> to SLDataSink */
void EnvelopeBunch::writeStats(FILE *f,double wfrac,OutputFileFormat oFormat) 
{
    static int counter  = 0, noHeader = 1, formatSDDS = 1;
    int nValid;
    double
        mc2e,
        ex,ey,
        b0,bRms,bMax,bMin,
        g0,dgdt,gInc,
        Rx,RxRms,RxMax,RxMin,
        Ry,RyRms,RyMax,RyMin,
        Px,PxRms,PxMax,PxMin,
        Py,PyRms,PyMax,PyMin,
        x0,x0Rms,x0Max,x0Min,
        y0,y0Rms,y0Max,y0Min,
        z0,zRms,zMax,zMin,
        I0,IRms,IMax,IMin;
    /*  Element
     *el;
     Field
     Ef,Bf; */

    if (noHeader) {
        switch (oFormat) {
        case oFormat_sdds:
            fprintf(f,"SDDS1\n");
            break;
        default :
            break;
        }

        fprintf(f,"&parameter name=wfrac, type=double, fixed_value=%20.14le, description=\"fractional window for statistics\" &end\n",
                wfrac);
        fprintf(f,"&parameter name=emtnx0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance x\" &end\n",1.0e6*emtnx0);
        fprintf(f,"&parameter name=emtny0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance y\" &end\n",1.0e6*emtny0);
        fprintf(f,"&parameter name=emtbx0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance Bz0x\" &end\n",1.0e6*emtbx0);
        fprintf(f,"&parameter name=emtby0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance Bz0y\" &end\n",1.0e6*emtby0);

        defineColumn(f,oFormat,"z",    "double","m",      "long. bunch position");
        defineColumn(f,oFormat,"t",    "double","ns",     "time");
        defineColumn(f,oFormat,"n",    "long",  NULL,     "number of used slices");
        defineColumn(f,oFormat,"E",    "double","MeV",    "beam energy");
        defineColumn(f,oFormat,"dEdt", "double","MeV/ps", "beam energy");
        defineColumn(f,oFormat,"dE",   "double","keV",    "rms energy spread");
        defineColumn(f,oFormat,"tau",  "double","ps",     "rms bunch length");
        defineColumn(f,oFormat,"Imax", "double","A",      "max bunch current");
        defineColumn(f,oFormat,"Irms", "double","A",      "rms bunch current");
        if (solver & sv_radial) {
            defineColumn(f,oFormat,"emtnx","double","mm mrad","normalized emittance x");
            defineColumn(f,oFormat,"emtny","double","mm mrad","normalized emittance y");
            defineColumn(f,oFormat,"Rx",   "double","mm",     "beam radius x");
            defineColumn(f,oFormat,"Ry",   "double","mm",     "beam radius y");
            defineColumn(f,oFormat,"RxMax","double","mm",     "Max. beam radius x");
            defineColumn(f,oFormat,"RyMax","double","mm",     "Max. beam radius y");
            defineColumn(f,oFormat,"RxMin","double","mm",     "Min. beam radius x");
            defineColumn(f,oFormat,"RyMin","double","mm",     "Min. beam radius y");
            defineColumn(f,oFormat,"Px",   "double","mrad",   "beam divergence x");
            defineColumn(f,oFormat,"Py",   "double","mrad",   "beam divergence y");
        }
        if (solver & sv_offaxis) {
            defineColumn(f,oFormat,"x0",   "double","mm",     "beam offset x");
            defineColumn(f,oFormat,"y0",   "double","mm",     "beam offset y");
            defineColumn(f,oFormat,"x0Max","double","mm",     "Max. offset x");
            defineColumn(f,oFormat,"y0Max","double","mm",     "Max. offset y");
            defineColumn(f,oFormat,"x0Min","double","mm",     "Min. offset x");
            defineColumn(f,oFormat,"y0Min","double","mm",     "Min. offset y");
        }
        defineColumn(f,oFormat,"dx",   "double","mm",     "offset of coordinate system x");
        defineColumn(f,oFormat,"dy",   "double","mm",     "offset of coordinate system y");
        defineColumn(f,oFormat,"dfi_x","double","deg",    "rotation of coordinate system x");
        defineColumn(f,oFormat,"dfi_y","double","deg",    "rotation of coordinate system y");
        defineColumn(f,oFormat,"Efield","double","MV/m",   "accelerating field");
        defineColumn(f,oFormat,"Bfield","double","T",      "field");

        switch (oFormat) {
        case oFormat_sdds :
            fprintf(f,"&data mode=ascii, no_row_counts=1,  &end\n");
            break;
        default :
            fprintf(f,"\n");
            break;
        }
        noHeader = 0;
    }

    runStats(sp_beta,wfrac,&b0,&bMax,&bMin,&bRms,&nValid);
    runStats(sp_z,wfrac,&z0,&zMax,&zMin,&zRms,&nValid);
    runStats(sp_I,wfrac,&I0,&IMax,&IMin,&IRms,&nValid);
    if (solver & sv_radial) {
        runStats(sp_Rx,wfrac,&Rx,&RxMax,&RxMin,&RxRms,&nValid);
        runStats(sp_Ry,wfrac,&Ry,&RyMax,&RyMin,&RyRms,&nValid);
        runStats(sp_Px,wfrac,&Px,&PxMax,&PxMin,&PxRms,&nValid);
        runStats(sp_Py,wfrac,&Py,&PyMax,&PyMin,&PyRms,&nValid);
        calcEmittance(wfrac,&ex,&ey,&nValid);
        //printf("emtny0 = %e/n/n",emtny0);
    }
    if (solver & sv_offaxis) {
        runStats(sp_x0,wfrac,&x0,&x0Max,&x0Min,&x0Rms,&nValid);
        runStats(sp_y0,wfrac,&y0,&y0Max,&y0Min,&y0Rms,&nValid);
    }
    calcEnergyChirp(wfrac,&g0,&dgdt,&gInc,&nValid);


    /** BET routine:  
     *
     * Ef = Field(); Bf = Field();
     * el = fe;
     * while (el) {
     *   Ef += el->getE(z0,t);
     *   Bf += el->getB(x0,y0,z0);
     *   el  = el->getNext();
     * }
     */

    /// Create dummies:

    double Efx=0,Efy=0,Efz=0,Bfx=0,Bfy=0,Bfz=0;

    /** calculate B and B field (average) for SDDS output
    */

    Bfz = AvBField();
    Efz = AvEField();

    mc2e = 1.0e-6*Physics::EMASS*Physics::c*Physics::c/Physics::q_e;

    switch (oFormat) {
    case oFormat_tab :
    case oFormat_sdds :
        // z, t
        fprintf(f,"%20.14le \t %20.14le", z0,1.0e9*t,z0);
        // n valid
        fprintf(f," \t %4d", nValid);
        // E, dEdt, Einc
        fprintf(f," \t %20.14le \t %20.14le \t %20.14le", mc2e*(g0-1.0),1.0e-12*mc2e*dgdt,1.0e3*mc2e*gInc);
        // tau
        fprintf(f," \t %20.14le", 1.0e12*zRms/Physics::c);
        // I, Irms
        fprintf(f," \t %20.14le \t %20.14le", IMax,Q*nValid*Physics::c/(zRms*sqrt(Physics::two_pi)*n));

        if (solver & sv_radial) {
            fprintf(f," \t %20.14le \t %20.14le",                  // emtnx, emtny
                    1.0e6*ex,1.0e6*ey);
            fprintf(f," \t %20.14le \t %20.14le",                  // Rx, Ry
                    1.0e3*Rx,1.0e3*Ry);
            fprintf(f," \t %20.14le \t %20.14le",                  // RxMax, RyMax
                    1.0e3*RxMax,1.0e3*RyMax);
            fprintf(f," \t %20.14le \t %20.14le",                  // RxMin, RyMin
                    1.0e3*RxMin,1.0e3*RyMin);
            fprintf(f," \t %20.14le \t %20.14le",                  // Px, Py
                    1.0e3*Px/Physics::c,1.0e3*Py/Physics::c);
        }

        if (solver & sv_offaxis) {
            fprintf(f," \t %20.14le \t %20.14le",                  // x0, y0
                    1.0e3*x0,1.0e3*y0);
            fprintf(f," \t %20.14le \t %20.14le",                  // x0Max, y0Max
                    1.0e3*x0Max,1.0e3*y0Max);
            fprintf(f," \t %20.14le \t %20.14le",                  // x0Min, y0Min
                    1.0e3*x0Min,1.0e3*y0Min);
        }

        // dx0, dy0
        fprintf(f," \t %20.14le \t %20.14le", 1.0e3*dx0, 1.0e3*dy0);
        // dfi_x, dfi_y
        fprintf(f," \t %20.14le \t %20.14le", 180.0*dfi_x/Physics::pi, 180.0*dfi_y/Physics::pi);
        // Ez, Bz
        fprintf(f," \t %20.14le \t %20.14le", 1.0e-6*(Efx+Efy+Efz), Bfx+Bfy+Bfz);
        fprintf(f,"\n");
        break;
    default :
        break;
    }
    fflush(f);
}

void EnvelopeBunch::writeSlice(FILE *f,OutputFileFormat oFormat) 
{
    Inform msg("writeSlice ");
    Inform msg2all("writeSlice ",INFORM_ALL_NODES);

    static int noHeader = 1, formatSDDS = 1;
    int i,nV;
    double CF,za,L;

    if (noHeader) {
        switch(oFormat) {
        case oFormat_sdds :
            fprintf(f,"SDDS1\n");
            fprintf(f,"&parameter name=time,   units=s,  type=double, description=\"Integration timestamp\" &end\n");
            fprintf(f,"&parameter name=z,      units=m,  type=double, description=\"Average position of bunch\" &end\n");
            fprintf(f,"&parameter name=dx0,    units=mm, type=double, description=\"Offset coordinate system x\" &end\n");
            fprintf(f,"&parameter name=dy0,    units=mm, type=double, description=\"Offset coordinate system y\" &end\n");
            fprintf(f,"&parameter name=dfix,   units=deg,type=double, description=\"Rotation coordinate system x\" &end\n");
            fprintf(f,"&parameter name=dfiy,   units=deg,type=double, description=\"Rotation coordinate system y\" &end\n");
            fprintf(f,"&parameter name=Q,      units=nC,type=double, fixed_value=%20.14le, description=\"EnvelopeBunch charge\" &end\n",1.0e9*Q);
            fprintf(f,"&parameter name=emtnx0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance x\" &end\n",1.0e6*emtnx0);
            fprintf(f,"&parameter name=emtny0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance y\" &end\n",1.0e6*emtny0);
            fprintf(f,"&parameter name=emtbx0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance Bz0x\" &end\n",1.0e6*emtbx0);
            fprintf(f,"&parameter name=emtby0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance Bz0y\" &end\n",1.0e6*emtby0);
            fprintf(f,"&parameter name=Bz0,    units=T,type=double, fixed_value=%20.14le, description=\"Magnetic field on cathode\" &end\n",Bz0);
            break;
        default :
            defineColumn(f,oFormat,"z","double","m","average position of bunch");
            break;
        }

        defineColumn(f,oFormat,"dz",   "double", "mm",   "slice position relative to center");
        defineColumn(f,oFormat,"t",    "double", "ps",   "slice position relative to center");
        defineColumn(f,oFormat,"E",    "double", "MeV",  "energy");
        if (solver & sv_radial) {
            defineColumn(f,oFormat,"Rx",   "double", "mm",   "slice radius x");
            defineColumn(f,oFormat,"Px",   "double", "mrad", "slice divergence");
        }
        if (solver & sv_offaxis) {
            defineColumn(f,oFormat,"x0",   "double", "mm",   "slice center x");
            defineColumn(f,oFormat,"Px0",  "double", "mrad", "transverse angle x");
        }
        if (solver & sv_radial) {
            defineColumn(f,oFormat,"Ry",   "double", "mm",   "slice radius y");
            defineColumn(f,oFormat,"Py",   "double", "mrad", "slice divergence");
        }
        if (solver & sv_offaxis) {
            defineColumn(f,oFormat,"y0",   "double", "mm",   "slice center y");
            defineColumn(f,oFormat,"Py0",  "double", "mrad", "transverse angle y");
        }
        defineColumn(f,oFormat,"Escz", "double", "MV/m", "slice longitudinal space-charge field");
        if (solver & sv_radial) {
            defineColumn(f,oFormat,"Escx", "double", "MV/m", "slice transverse space-charge field x");
            defineColumn(f,oFormat,"Escy", "double", "MV/m", "slice transverse space-charge field y");
        }
        if (solver & sv_offaxis) {
            defineColumn(f,oFormat,"Exw",  "double", "MV/m", "slice transverse wake-field x");
            defineColumn(f,oFormat,"Eyw",  "double", "MV/m", "slice transverse wake-field y");
        }
        defineColumn(f,oFormat,"Ezw",  "double", "MV/m", "slice longitudinal wake-field");
        defineColumn(f,oFormat,"I",    "double", "A",    "beam current");
        if (solver & sv_offaxis) {
            defineColumn(f,oFormat,"Ex",   "double", "MV/m", "Eacc field x");
            defineColumn(f,oFormat,"Ey",   "double", "MV/m", "Eacc field y");
        }
        defineColumn(f,oFormat,"Ez",   "double", "MV/m", "Eacc field z");
        defineColumn(f,oFormat,"Bx",   "double", "T",    "B field x");
        defineColumn(f,oFormat,"By",   "double", "T",    "B field y");
        defineColumn(f,oFormat,"Bz",   "double", "T",    "B field z");
        if (solver & sv_radial) {
            defineColumn(f,oFormat,"Kx",   "double", "1/m2", "Focus strength x");
            defineColumn(f,oFormat,"Ky",   "double", "1/m2", "Focus strength y");
            defineColumn(f,oFormat,"Kz",   "double", "1/m2", "Focus strength z");
        }
        if (solver & sv_offaxis) {
            defineColumn(f,oFormat,"Dx",   "double", "1/m",  "Deflection x");
            defineColumn(f,oFormat,"Dy",   "double", "1/m",  "Deflection y");
            defineColumn(f,oFormat,"Dz",   "double", "1/m",  "Deflection z");
        }

        switch (oFormat) {
        case oFormat_sdds :
            fprintf(f,"&data mode=ascii &end\n");
            break;
        default:
            fprintf(f,"\n");
            break;
        }
        noHeader = 0;
    }

    za = zAvg();
    L  = zHead() - zTail();
    CF = 0.0;

    if (s) { /// there are slices defined
        CF = (Q==0.0?0.0:0.25*(I->max()/Physics::Ia)*Physics::EMASS*Physics::c*Physics::c/Physics::q_e);
    }

    nV = 0;
    for (i=0; i<n; i++) {
        if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) ++nV;
    }
    if (nV > 0) {
        switch (oFormat) {
        case oFormat_sdds :    
            fprintf(f,"! next page\n");
            fprintf(f,"%20.12le\n",t);
            fprintf(f,"%20.12le\n",za);
            fprintf(f,"%20.12le\n",1.0e3*dx0);
            fprintf(f,"%20.12le\n",1.0e3*dy0);
            fprintf(f,"%20.12le\n",180.0*dfi_x/Physics::pi);
            fprintf(f,"%20.12le\n",180.0*dfi_y/Physics::pi);
            fprintf(f,"%d\n",nV);
        case oFormat_tab :
            for (i=0; i<n; i++) {
                double c2 = Physics::c*Physics::c;

                if ((s[i]->p[SLI_z] > zCat) && s[i]->valid) {
                    double
                        *v    = s[i]->p,
                        gamma = s[i]->gamma();
                    //	  Element *el;

                    if (oFormat != oFormat_sdds) {
                        fprintf(f,"%20.12le \t",za);
                    }
                    //dz, t
                    fprintf(f,"%20.12le \t %20.12le ", 1.0e+03*(v[SLI_z] - za),1.0e+12*(v[SLI_z] - za)/(Physics::c*v[SLI_beta]));
                    //E, dE
                    fprintf(f,"\t %20.12le ",0.511*(gamma-1.0));
                    if (solver & sv_radial) {
                        //Rx, Px
                        fprintf(f,"\t %20.12le \t %20.12le ", 2.0e+03*v[SLI_x],2.0e+03*v[SLI_px]/(v[SLI_beta]*Physics::c));
                    }
                    if (solver & sv_offaxis) {
                        //x0, px0
                        fprintf(f,"\t %20.12le \t %20.12le ", 1.0e+03*v[SLI_x0],1.0e+03*v[SLI_px0]/Physics::c);
                    }
                    if (solver & sv_radial) {
                        //Ry, Py
                        fprintf(f,"\t %20.12le \t %20.12le ", 2.0e+03*v[SLI_y],2.0e+03*v[SLI_py]/(v[SLI_beta]*Physics::c));
                    }
                    if (solver & sv_offaxis) {
                        //y0, py0
                        fprintf(f,"\t %20.12le \t %20.12le ", 1.0e+03*v[SLI_y0],1.0e+03*v[SLI_py0]/Physics::c);
                    }

                    //Esct
                    fprintf(f,"\t %20.12le ", 1.0e-06*Esct[i]);
                    if (solver & sv_radial) {
                        //Escx Escy
                        fprintf(f,"\t %20.12le \t %20.12le ", 1.0e-06*CF*G[i]/(v[SLI_x]*v[SLI_beta]), 1.0e-06*CF*G[i]/(v[SLI_y]*v[SLI_beta]));
                    }
                    if (solver & sv_offaxis) {
                        //Exw, Eyw
                        fprintf(f,"\t %20.12le \t %20.12le ", 1.0e-6*Exw[i],1.0e-6*Eyw[i]);
                    }

                    //Ezw
                    fprintf(f,"\t %20.12le ", 1.0e-6*Ezw[i]);

                    // I
                    fprintf(f,"\t %20.12le ", Q>0.0?I->get(v[SLI_z],itype_lin):0.0);

                    /** access field data by BF and EF:
                     *
                     * BF[i] -> gives one Vector_t with three components from the
                     *      ith slice
                     *
                     * BF[i](j) -> accesses the jth component of the ith slice
                     *
                     * Same with EF. 
                     */


                    if (solver & sv_offaxis) {
                        // Ex, Ey
                        fprintf(f,"\t %20.12le \t %20.12le ", 1.0e-6*EF[i](0),1.0e-6*EF[i](1));
                    }

                    // Ez
                    fprintf(f,"\t %20.12le ", 1.0e-6*EF[i](2));
                    // Bx, By, Bz
                    fprintf(f,"\t %20.12le \t %20.12le \t %20.12le ", BF[i](0),BF[i](1),BF[i](2));

                    if (solver & sv_radial) {
                        fprintf(f,"\t %20.12le \t %20.12le \t %20.12le ",   // Kx, Ky, Kz
                                KR[i](0)/c2,KR[i](1)/c2,KR[i](2)/c2);
                    }
                    if (solver & sv_offaxis) {
                        fprintf(f,"\t %20.12le \t %20.12le \t %20.12le ",   // Dx, Dy, Dz
                                KT[i](0)/c2,KT[i](1)/c2,KT[i](2)/c2);
                    }
                    fprintf(f,"\n");
                }
            }
            break;
cefault :
            break;
        }
        //msg << "finished function writeSlice" << endl;
    }
}



/// z [m] of center 
/// width (see .h file)
/// Rect: slope of rect, Gauss: fraction used
void EnvelopeBunch::setLShape(EnvelopeBunchShape shape, double z0, double w, double frac)   
{
    int n2 = n/2, i;
    double sqr2 = sqrt(2.0),dz0,sl,v,vMax;

    switch (shape) {
    case bsRect :
        sl   = fabs(frac) > fabs(w) ? w/2.0 : frac;
        double ffact = w < 0.0 ? Physics::c : 1.0;
        vMax = ffact*(fabs(w) - fabs(sl))/2.0;
        dz0  = ffact*(fabs(w) + 0.95*fabs(sl))/(n-1);
        sl   = ffact*fabs(sl);
        v    = 0.0;
        s[n2]->p[SLI_z] = z0;

        for (i=1; i<=n/2; i++) {
            if (v <= vMax) {
                v += dz0;
            } else {
                v += (2.0*dz0/(1.0+cos(Physics::u_two_pi*(v-vMax)/sl)));
            }
            if (n2+i <  n) 
                s[n2+i]->p[SLI_z] = z0 + v*s[n2+i]->p[SLI_beta];
            if (n2-i >= 0) 
                s[n2-i]->p[SLI_z] = z0 - v*s[n2-i]->p[SLI_beta];

            //      s[i]->p[SLI_z] = z0 + (fabs(w)*((1.0*i/(n-1)) - 0.5)*(w<0.0?(c*s[i]->p[SLI_beta]):1.0));
        }

        I0avg = Q*(w<0.0?1.0:Physics::c)/fabs(2.0*w);
        break;
    case bsGauss :
        s[n2]->p[SLI_z] = z0;
        for (i=1; i<=n/2; i++) {
            rootValue = 1.0 - 2.0*i*frac/(n+1);
            v         = fabs(w)*sqr2*findRoot(erfRoot,0.0,5.0,1.0e-5)*(w<0.0?Physics::c:1.0);
            if (n2+i <  n) 
                s[n2+i]->p[SLI_z] = z0 + v*s[n2+i]->p[SLI_beta];
            if (n2-i >= 0) 
                s[n2-i]->p[SLI_z] = z0 - v*s[n2-i]->p[SLI_beta];
        }
        I0avg = 0.0;
        break;
    }

    backup();  /// save new state
}

 /** normalized emittance (x,y) [m rad]
 */
 /** transverse size (x,y) [m]
 */
void EnvelopeBunch::setTShape(double enx, double eny, double rx, double ry, double b0) 
{
    int i;

    /// set emittances
    emtnx0 = enx; 
    emtny0 = eny;
    emtbx0 = Physics::q_e*rx*rx*Bz0/(8.0*Physics::EMASS*Physics::c);
    emtby0 = Physics::q_e*ry*ry*Bz0/(8.0*Physics::EMASS*Physics::c);

    /// set the diameter of each slice
    for (i=0; i<n; i++) {
        s[i]->p[SLI_x] = rx/2.0;
        s[i]->p[SLI_px] = 0.0;
        s[i]->p[SLI_y] = ry/2.0;
        s[i]->p[SLI_py] = 0.0;
    }

    /// save new state
    backup();
}

/// horizontal centroid + angle
/// vertical centroid + angle
void EnvelopeBunch::setTOffset(double x0, double px0, double y0, double py0) 
{
    for(int i=0; i<n; i++) {
        s[i]->p[SLI_x0] = x0; s[i]->p[SLI_px0] = px0;
        s[i]->p[SLI_y0] = y0; s[i]->p[SLI_py0] = py0;
    }
}

 /// mean energy (eV)
 /// correlated energy spread [eV/m]
void EnvelopeBunch::setEnergy(double E0, double dE) 
{
    double g0    = 1.0 + (Physics::q_e*E0/(Physics::EMASS*Physics::c*Physics::c));
    double dg    = fabs(dE)*Physics::q_e/(Physics::EMASS*Physics::c*Physics::c);
    double z0    = zAvg();
    for(int i=0; i<n; i++) {
        double g = g0 + (s[i]->p[SLI_z] - z0)*dg;
        s[i]->p[SLI_beta] = sqrt(1.0 - (1.0/(g*g)));
    }

    /// save new states
    backup();
}

double EnvelopeBunch::getGamma(int i) 
{
    //double tempG;
    //tempG = 1/sqrt(1-(s[i]->p[SLI_beta])*(s[i]->p[SLI_beta]));
    //return tempG;
    return s[i]->gamma();
}

double EnvelopeBunch::getBeta(int i) 
{
    return s[i]->p[SLI_beta];
}

/** Set emittances 
*/

void EnvelopeBunch::setEx(double emi) 
{
    emtnx0=emi;
}

void EnvelopeBunch::setEy(double emi) 
{
    emtny0=emi;
}

/** Output of Rx 
*/

void EnvelopeBunch::plotR() 
{
    for(int i=0;i<getN();i++) {
        printf("rx[%d] = %e\n",i,s[i]->p[SLI_x]);
    }
}

/** Calculates the average (over all slices) Bfield from BF 
*/

double EnvelopeBunch::AvBField() 
{
    double bf=0.0;
    for(int slice=0; slice<getN(); slice++) {
        for(int i=0; i<3; i++) {
            bf += BF[slice](i);
        }
    }
    return bf/getN();
}

/** Calculates the average (over all slices) Efield from EF 
*/

double EnvelopeBunch::AvEField() 
{
    double ef=0.0;
    for(int slice=0; slice<getN(); slice++) {
        for(int i=0; i<3; i++) {
            ef += EF[slice](i);
        }
    }
    return ef/getN();
}

/** sets the charge. "sign" gives the sign of Q, since BET only accepts 
 *  positive values for Q 
 */

void EnvelopeBunch::setCharge(double _Q) 
{
    sign = _Q<0.0?-1:1;
    Q = abs(_Q);
}

/** \brief set solver 
 * just sets the solver, all error messages were deleted
 */

void EnvelopeBunch::setSolver(int s) 
{
    solver = s;
    /**writeError(errModeMaster,errMessage,"Differential Equation (DE) solving method: %s",
      solver&sv_fixedStep?"Fixed step 4th order RK":"Variable step RK");
      writeError(errModeMaster,errMessage,"Field evaluation: %s the DE solver",
      solver&sv_fieldOutside?"Performed before calling":"Inside");
      writeError(errModeMaster,errMessage,"Transverse dynamics: %s in DV",
      solver&sv_radial?"included":"EXCLUDED");
      writeError(errModeMaster,errMessage,"Off-axis dynamics:   %s in DV",
      solver&sv_offaxis?"included":"EXCLUDED");
      writeError(errModeMaster,errMessage,"Longitudinal wakes:  %s in DV",
      solver&sv_lwakes?"included":"EXCLUDED");
      writeError(errModeMaster,errMessage,"Transverse wakes:    %s in DV",
      solver&sv_twakes?"included":"EXCLUDED");
      writeError(errModeMaster,errMessage,"Path tracking:       %s",
      solver&sv_s_path?"s-axis":"z-axis");
      */
}


/** \brief cWake()
 * Calculate wake fields over whole bunch, this function got deleted
 */

/* 
void EnvelopeBunch::cWake() 
{
   int 
   i,i0,i_step,j,nV;
   double
   qc,z0,g0;
   Element 
 *ep;
 Slice
 *cs;

 if (n < 2) {
 writeError(errModeMaster,errGeneral,
 "EnvelopeBunch::cWake() called with insufficient slices (%d)",n);
 return;
 }

#ifdef USE_MPI
i0 = mpi_rank; i_step = mpi_size;
#else 
i0 = 0; i_step = 1;
#endif

for (i=0; i<n; ++i) {
Ezw[i] = 0.0;
Exw[i] = 0.0;
Eyw[i] = 0.0;
}

nV = 0;

if ((Q > 0.0) && (solver & (sv_lwakes | sv_twakes))) {
g0 = 0.0;
for (i=0; i<n; ++i) {
cs = s[i];
if (cs->p[SLI_z] > zCat) {
g0 += cs->gamma();
++nV;
}
}
}

if (nV > 1) {
g0 = g0/nV;
qc = Q/(c*n);
for (i=i0; i<n; i += i_step) { // loop over slices
cs = s[i];
if (cs->p[SLI_z] > zCat) {
ep = fE;
z0  = cs->p[SLI_z];
while (ep) { // loop over RF structures
switch (ep->getType()) {
case etRF:
case etSWA:
case etTWA:
if (ep->isActive(z0) && 
(((RF *) ep)->getIris() > 0.0)) {
for (j=0; j<n; j++) {
if (s[j]->p[SLI_z] > z0) {
Slice  *sp;
double wt,g;

sp = s[j];
g  = sp->gamma();
wt = (g - g0)*((RF *)ep)->getWakeT(sp->p[SLI_z] - z0)/g;
if (solver & sv_twakes) {
Exw[i] -= (wt*(sp->p[SLI_x0] - ep->get(ets_x0)));
Eyw[i] -= (wt*(sp->p[SLI_y0] - ep->get(ets_y0)));
}
if (solver & sv_lwakes) {
Ezw[i] -= (((RF *)ep)->getWakeL(sp->p[SLI_z] - z0));
}
}
}
Exw[i] *= qc;
Eyw[i] *= qc;
Ezw[i] *= qc;
}
break;
default:
break;
}
ep = ep->getNext();
}
} 
} 
} 
} 




EnvelopeBunch::calcI()
    calculates the current profile irrespective of a reshuffling of the
    slice position.
    If the bunch charge is <= 0.0, a bunch charge of e is assumed



#ifdef ALT_CALCI
    */ /*
          void EnvelopeBunch::calcI() {
          static int
          already_called = 0;

          Profile
        *qp;
        double
        z0,dt2,dt_min,
        tMin,tMax,dt,
        *tau,*t1,*t2,*q,*tq,*z1,*I1;
        int
        i,j,Nb;

        if ((dStat & ds_currentCalculated) || 
        (already_called && (Q <= 0.0)) ) return;

        already_called = 1;

    // delete the old profile
    if (I) delete I;

    /* step 1:
    calculate the average bunch postion

    z0 = 0.0;
    for (i=0; i<n; i++) {
    z0 += s[i]->p[SLI_z];
    }
    z0 = z0/n;

    /* step 2:
    convert the position of a slice to time (tau)

    tau = vector(n);
    for (i=0; i<n; i++) {
    tau[i] = (s[i]->p[SLI_z] - z0)/(s[i]->p[SLI_beta]*c);
    }

    /* step 3:
    - calculate the start (t1) and stop (t2) time of each slice
    - calculate the average slice duration

    t1  = vector(n); t2 = vector(n);
    dt2 = 0.0;
    t1[0] = 0.5*(tau[0] - tau[1]);
    t2[0] = 0.5*(tau[0] + tau[1]);
    dt2   = pow(t2[0]-t1[0],2);
    for (i=1; i<n-1; i++) {
    t1[i] = t2[i-1];
    t2[i] = 0.5*(tau[i+1] + tau[i]);
    dt2  += pow(t2[i]-t1[i],2);
    }
    t1[n-1] = t2[n-2];
    t2[n-1] = 0.5*(3.0*tau[n-1] - tau[n-2]);
    dt2    += pow(t2[n-1]-t1[n-1],2);

    /* step 3:
    - make sure that t1[i] < t2[i]
    - get bunch duration (tMin,tMax)
    - spike suppression 1:
    increase the slice duration if it happens to be too short:
    heuristic approach: dt > dt_avg/5 

    dt_min = 0.2*sqrt(dt2/n);
    for (i=0; i<n; i++) {
    if (t2[i] < t1[i]) {
    swap(&t1[i],&t2[i]);
    }
    if ((t2[i]-t1[i]) < dt_min) {
    t1[i] = tau[i] - dt_min/2.0;
    t2[i] = tau[i] + dt_min/2.0;
    }
    if (i == 0) {
    tMin = t1[0]; tMax = t2[0];
    } else {
    if (t1[i] < tMin) tMin = t1[i];
    if (t2[i] > tMax) tMax = t2[i];
    }
    }

    step 4:
    - bin charge profile in q (Nb bins)
    - transform the histogram into a profile 

    Nb = (int) ((tMax - tMin)/dt_min);
    if (Nb < 100)   Nb = 100;
    q  = vector(Nb);
    tq = vector(Nb);
    dt = (tMax - tMin)/(Nb - 1);
    for (j=0; j<Nb; j++) {
    q[j]  = 0.0;
    tq[j] = tMin + dt*(2*j + 1)/2.0;
    }
    for (i=0; i<n; i++) {
    int
    j1,j2;

    j1 = (int) ((t1[i] - tMin)/dt);
    j2 = (int) ((t2[i] - tMin)/dt);
    if (j1 == j2) {
    q[j1] += (Q/n);
    } else {
    double
    f1,f2;

    // fraction in first and last bin
    f1 = (tMin + dt*(j1+1) - t1[i])/(t2[i]-t1[i]);
    f2 = (t2[i] - (tMin + dt*(j2)))/(t2[i]-t1[i]);
    q[j1] += (f1*Q/n);
    q[j2] += (f2*Q/n);
    if ((j2-j1) > 1) {
    double 
    dq;

    dq = (1.0 - f1 - f2)*Q/n/(j2-j1-1);
    for (j=j1+1; j<j2; j++) q[j] += dq;
    }
    }
    }
    sgSmooth(q,Nb,10,10,0,1);
    qp = new Profile(tq,q,Nb);

    /* step 5:
    transform the charge profile into a current profile 

    z1 = vector(n);
    I1 = vector(n);
    for (i=0; i<n; i++) {
    z1[i] = s[i]->p[SLI_z];
    I1[i] = qp->get(tau[i])/dt; 
    }
    I = new Profile(z1,I1,n);

    free(tau); free(t1); free(t2);
    free(tq);  free(q);
    free(z1);  free(I1);
    delete qp;
    } // calcI() alternative method
    */
    /*
#else
*/

    /** \brief calculates current 
     *  function to calculate the current profile, the second of two programmes 
     */

void EnvelopeBunch::calcI() 
{

    Inform msg("calcI ");
    double temp;

    static int already_called = 0;

    int i,j,k,np,n1,n2;
    double beta,bSum,Mz1,MI1,dz,dzMin,dz2Sum,sigma_dz,z,zMin,zMax,q,Qcalc,*b,*z1,*z2,*I1,*I2;

    if ((dStat & ds_currentCalculated) || 
            (already_called && (Q <= 0.0)) ) return;

    already_called = 1;

    // delete the old profile
    if (I) delete I;

    z1     = vector(n);
    b      = vector(n);
    bSum   = 0.0;
    dz2Sum = 0.0;
    n1     = 0;
    for (i=0; i<n; i++) {
        beta   = s[i]->p[SLI_beta];
        if (beta > 0.0) {
            b[n1]   = beta;
            z1[n1]  = s[i]->p[SLI_z];
            if (n1 > 0) dz2Sum += pow(z1[n1]-z1[n1-1],2);
            bSum  += beta;
            ++n1;
        } 
    }
    if (n1 < 2) {
        free(z1);
        free(b);
        //writeError(errModeMaster,errWarning,
        //       "Insufficient points to calculate the current (n1 = %d)",
        //       n1);
        writeError();
        I = new Profile(0.0);
        return;
    }

    sigma_dz = sqrt(dz2Sum/(n1-1));
    beta = bSum/n1;
    sort2(z1,b,n1);

    q = Q>0.0?Q/n:Physics::q_e;

    /// 1. determine current from the slice distance

    I1 = vector(n1);
    I1[0]  = I1[n1-1] = 0.0;

    ///limit the max current to 5x the sigma value
    /// to reduce noise problems 
    dzMin  = 0.2*sigma_dz;
    for (i=1; i<n1-1; i++) {
        j = 1;
        do {
            dz = fabs(z1[i+j]-z1[i-j]);
            ++j;
        } while ((dz < dzMin*(j-1)) && ((i+j) < n1) && ((i-j) >= 0));

        if ((dz >= dzMin*(j-1)) && ((i+j) < n1) && ((i-j) >= 0)) {
            I1[i] = 0.25*q*Physics::c*(b[i+j]+b[i-j])/(dz*(j-1));
        } else {
            I1[i] = I1[i-1];
        }
    }
    I1[0]    = I1[1];
    I1[n1-1] = I1[n1-2]; 

    ///2. Remove points with identical z-value
    /// and then smooth the current profile

    zMin = zTail();
    zMax = zHead();
    dz = (zMax - zMin)/n; // create a window of the average slice distance
    z2 = vector(n1);
    I2 = vector(n1);
    Mz1= 0.0;
    MI1= 0.0;
    np = 0;

    /// first value
    j = 0;
    while ((j<n1) && ((z1[j]-z1[0]) <= dz)) {
        Mz1 += z1[j];
        MI1 += I1[j];
        ++j;
        ++np;
    }
    z2[0] = Mz1/np;
    I2[0] = MI1/np;

    /// following values
    k = 0;
    for (i=1; i<n1; i++) {
        /// add new points
        j = 0;
        while (((i+j)<n1) && ((z1[i+j]-z1[i]) <= dz)) {
            if ((z1[i+j]-z1[i-1]) > dz) {
                Mz1 += z1[i+j];
                MI1 += I1[i+j]; 
                ++np;
            }
            ++j;
        }

        // remove obsolete points
        j = 1;
        while (((i-j)>=0) && ((z1[i-1]-z1[i-j]) < dz)) {
            if ((z1[i]-z1[i-j]) > dz) {
                Mz1 -= z1[i-j];
                MI1 -= I1[i-j]; 
                --np;
            }
            ++j;
        } 
        z2[i-k] = Mz1/np;
        I2[i-k] = MI1/np;

        /// make sure there are no duplicate z coordinates
        if (z2[i-k] <= z2[i-k-1]) {
            I2[i-k-1] = 0.5*(I2[i-k] + I2[i-k-1]);
            ++k;
        }
    }
    n2 = n1 - k;
    if (n2 < 1) {
        *gmsg << "Insufficient points to calculate the current (m = " << n2 << ")" << endl;
        I = new Profile(0.0);

    } else {
        /// 3. smooth further
        if (n2 > 40) {
            sgSmooth(I2,n2,n2/20,n2/20,0,1);
        }

        /// 4. create current profile 
        I = new Profile(z2,I2,n2);

        /**5. Normalize profile to match bunch charge as a constant
         * However, only normalize for sufficient beam energy
         */

        itsBunch = this; 

        Qcalc    = 0.0; 
        z   = zMin;
        dz  = (zMax - zMin)/99;
        for (i=1; i<100; i++) {
            Qcalc += I->get(z,itype_lin);
            z += dz;
        }
        Qcalc *= (dz/(beta*Physics::c));
        I->scale((Q>0.0?Q:Physics::q_e)/Qcalc);
        temp = Qcalc;
    }
    free(z2); free(I2);

    free(b);
    free(z1); free(I1);

    dStat |= ds_currentCalculated;
}

//#endif

/** \brief cSpaceCharge()
 * Calculate longitudinal-, and transverse space-charge fields for bunch
 * Output is stored in the global class arrays Esct and G, respectively.
 */

void EnvelopeBunch::cSpaceCharge() 
{
    string fn("SCValues_OPAL.dat");
    Inform msg ("cSpaceCharge");
    int i,i0,i_step,j,nV;
    double A0,Aj,Imax,dz,dzMin,z0,zj,v,sm,*xi,*z;
    EnvelopeSlice *cs;

    if (n < 2) {
        *gmsg << "EnvelopeBunch::cSpaceCharge called with insufficient slices (" << n << ")" << endl;
        return;
    }

//#ifdef USE_MPI
    //i0 = mpi_rank; 
    //i_step = mpi_size;
//#else
    //i0 = 0; 
    //i_step = 1;
//#endif 

    i0=0;

    for (i=0; i<n; ++i) {
        Esct[i] = 0.0;
        G[i]    = 0.0;
    }

    if ((Q <= 0.0) || (I->max() <= 0.0)) {
        msg << "going out" << endl;
        return;
    }

    /// 1. prepare variables

    Imax  = I->max();
    xi    = vector(n);
    z     = vector(n);
    nV    = 0;
    sm    = 0.0;

    for (i=0; i<n; ++i) {
        cs = s[i];
        z[i]  = cs->p[SLI_z];
        if (cs->p[SLI_z] > zCat) {
            ++nV;
            A0    = 4.0*cs->p[SLI_x]*cs->p[SLI_y];
            sm   += A0;
            xi[i] = A0*(1.0 - pow(cs->p[SLI_beta],2)); // g2
        }
    }

    if (nV < 2) {
        free(xi); free(z);
        return;
    }
    A0 = sm/nV;

    /// 2. calculate the space charge

    dzMin = 5.0*Physics::c*Q/(Imax*n);
    for (i=i0; i<n; i++) {

        /** Only include space-charge for particles 
         * that are downstream of the cathode 
         */

        //if ((z[i] > zCat) && (s[i]->p[SLI_beta] > BETA_MIN1)) {
        double bt;

        z0 = z[i];

        sm = 0.0;
        for (j=0; j<n; j++) {
            cs = s[j];
            zj = cs->p[SLI_z];
            dz = fabs(zj - z0);
            if ((dz > dzMin) && (zj > zCat)) {
                Aj = xi[j]/pow(dz,2);
                v  = 1.0 - (1.0/sqrt(1.0 + Aj));
                if (zj > z0) {
                    sm -= v;
                } else {
                    sm += v;
                }
            }
        }

        /// longitudinal effect

        bt = s[i]->p[SLI_beta];
        Esct[i] = (bt<BETA_MIN1?0.0:bt<BETA_MIN2?pow((bt-BETA_MIN1)/(BETA_MIN2-BETA_MIN1),2):1.0)*Q*sm/(Physics::two_pi*Physics::epsilon_0*A0*(nV-1)); 
        G[i] = I->get(z0,itype_lin)/Imax;
        if (bt < BETA_MIN2) { // tweak to compensate for non-relativity
            if (s[i]->p[SLI_beta] < BETA_MIN1) {
                G[i] = 0.0;
            } else {
                G[i] *= pow((bt-BETA_MIN1)/(BETA_MIN2-BETA_MIN1),2);
            }
        }
        printSC(fn.c_str());
    }
    free(z);
    free(xi);
}

//ff
/** function to write spacecharge (SC) output into file 
*/

void EnvelopeBunch::printSC(const char *fName) 
{

    FILE* fout = NULL;

    fout = fopen(fName,"w");
    fprintf(fout,"#Output for space charge \n#i\tEsct[i]\tG[i]\n");
    for(int i=0;i<10;i++) {
        fprintf(fout,"%d\t%e\t%e\n",i,Esct[i],G[i]);
    }
    fclose(fout);
}

//ff
void EnvelopeBunch::writeSC(int slices) 
{

    Inform msg ("writeSC");

    msg << "Slice no i\tG[i]\tEsct" << endl;
    for(int i=0;i<slices;i++) {
        msg << i << "\t" << G[i] << "\t" << Esct[i] << endl;
    }  
}


double EnvelopeBunch::moveZ0(double zC) 
{
    double dz;
    int    i;

    zCat = zC;
    dz   = zC - zHead();
    if(dz > 0.0) {
        for(i=0; i<n; i++) {
            s[i]->p[SLI_z] += dz;
        }
        backup(); // save the new state 
        //writeError(errModeMaster,errMessage,
        //       "EnvelopeBunch::moveZ0(): bunch moved with %.3f m to %.3f m",
        //       dz,zCat);
        writeError();
    }

    return dz;
}

void EnvelopeBunch::backup() 
{
    for(int i = 0; i<n; i++) 
        s[i]->backup();
}


/** EnvelopeBunch::checkScreen check if present position is marked as a screen
 * returns 1 if all slices have passed through the screen 
 *
 *  Note: actual dumping of screen may be delayed if space-charge fields
 *        are not updated on each time-step. THIS ROUTINE ONLY BECOMES
 *        ACTIVE AFTER A SYNCHRONIZATION OF THE SLICES, WHICH ONLY HAPPENS
 *        AFTER A SPECIFIED NUMBER OF TIME-STEPS (DEFAULT=1) OR AFTER A 
 *        TIME OR QUICK OUTPUT. SEE CODE of system::run() for details
 */
int EnvelopeBunch::checkScreen(double zScreen,FILE *f, OutputFileFormat oFormat) 
{
    static int noHeader = 1;

    int i,nV,allMarked = 1;

    nV = 0;
    for(i=0; i<n; i++) {
        if(s[i]->valid) {
            int result;

            /// I have to mark all before parsing he error!
            result = s[i]->mark(zScreen,t);
            if(allMarked == 1) 
                allMarked = result;
            ++nV;
        }
    }
    if((allMarked == 1) && (nV > 0)) {
        double t0;

        if(noHeader) {
            switch(oFormat) {
            case oFormat_sdds :
                fprintf(f,"SDDS1\n");
                fprintf(f,"&parameter name=time,   units=s, type=double, description=\"average arrival time\" &end\n");
                fprintf(f,"&parameter name=z,      units=m, type=double, description=\"Screen position\" &end\n");
                fprintf(f,"&parameter name=Q,      units=nC,type=double, fixed_value=%20.14le, description=\"EnvelopeBunch charge\" &end\n",1.0e9*Q);
                fprintf(f,"&parameter name=emtnx0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance x\" &end\n",1.0e6*emtnx0);
                fprintf(f,"&parameter name=emtny0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance y\" &end\n",1.0e6*emtny0);
                fprintf(f,"&parameter name=emtbx0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance Bz0x\" &end\n",1.0e6*emtbx0);
                fprintf(f,"&parameter name=emtby0, units=\"mm mrad\",type=double, fixed_value=%20.14le, description=\"Norm. Intrinsic Emittance Bz0y\" &end\n",1.0e6*emtby0);
                fprintf(f,"&parameter name=Bz0,    units=T,type=double, fixed_value=%20.14le, description=\"Magnetic field on cathode\" &end\n",Bz0);
                break;
            default :
                defineColumn(f,oFormat,"z","double","m","screen position");
                break;
            }

            defineColumn(f,oFormat,"t",    "double", "ps",   "delay-time of arrival");
            defineColumn(f,oFormat,"E",    "double", "MeV",  "energy");
            if (solver & sv_radial) {
                defineColumn(f,oFormat,"Rx",   "double", "mm",   "slice radius x");
                defineColumn(f,oFormat,"Px",   "double", "mrad", "slice divergence");
            }
            if (solver & sv_offaxis) {
                defineColumn(f,oFormat,"x0",   "double", "mm",   "slice center x");
                defineColumn(f,oFormat,"Px0",  "double", "mrad", "transverse angle x");
            }
            if (solver & sv_radial) {
                defineColumn(f,oFormat,"Ry",   "double", "mm",   "slice radius y");
                defineColumn(f,oFormat,"Py",   "double", "mrad", "slice divergence");
            }
            if (solver & sv_offaxis) {
                defineColumn(f,oFormat,"y0",   "double", "mm",   "slice center y");
                defineColumn(f,oFormat,"Py0",  "double", "mrad", "transverse angle y");
            }

            switch (oFormat) {
            case oFormat_sdds :
                fprintf(f,"&data mode=ascii &end\n");
                break;
            default:
                fprintf(f,"\n");
                break;
            }
            noHeader = 0;
        }

        t0 = 0.0;
        for (i=0; i<n; i++) {
            if (s[i]->valid) {
                t0 += (s[i]->p_scr[SLI_z]/nV);
            }
        }

        switch (oFormat) {
        case oFormat_sdds :
            fprintf(f,"! next page\n");
            fprintf(f,"%20.12le\n",t0);
            fprintf(f,"%20.12le\n",zScreen);
            fprintf(f,"%d\n",nV);
        case oFormat_tab :
            for (i=0; i<n; i++) {
                if (s[i]->valid) {
                    double 
                        *v    = s[i]->p_scr, 
                        gamma = 1.0/sqrt(1.0 - pow(v[SLI_beta],2));

                    if (oFormat != oFormat_sdds) {
                        fprintf(f,"%20.12le \t",zScreen);
                    }

                    fprintf(f,"\t %20.12le ",                             // dt (negative -> like z-screens)
                            1.0e+12*(t0 - v[SLI_z]));
                    fprintf(f,"\t %20.12le ",0.511*(gamma-1.0));          //E, dE
                    if (solver & sv_radial) {
                        fprintf(f,"\t %20.12le \t %20.12le ",               //Rx, Px
                                2.0e+03*v[SLI_x],2.0e+03*v[SLI_px]/(v[SLI_beta]*Physics::c));
                    }
                    if (solver & sv_offaxis) {
                        fprintf(f,"\t %20.12le \t %20.12le ",               //x0, px0
                                1.0e+03*v[SLI_x0],1.0e+03*v[SLI_px0]/(v[SLI_beta]*Physics::c));
                    }
                    if (solver & sv_radial) {
                        fprintf(f,"\t %20.12le \t %20.12le ",               //Ry, Py
                                2.0e+03*v[SLI_y],2.0e+03*v[SLI_py]/(v[SLI_beta]*Physics::c));
                    }
                    if (solver & sv_offaxis) {
                        fprintf(f,"\t %20.12le \t %20.12le ",               //y0, py0
                                1.0e+03*v[SLI_y0],1.0e+03*v[SLI_py0]/(v[SLI_beta]*Physics::c));
                    }
                    fprintf(f,"\n");
                }
            }
            break;
cefault :
            break;
        }
    }

    if (allMarked != 0) { // clear screen for next use
        for (i=0; i<n; i++) {
            if (s[i]->p_scr) free(s[i]->p_scr);
            s[i]->p_scr = NULL;
        }
    }

    return allMarked;
}

/** tReset()
 *  reset the timer of the bunch to zero 
 */
double EnvelopeBunch::tReset(double dt) 
{
    double new_dt = dt;

    if(dt == 0.0) {
        new_dt = t;
        //writeError(errModeMaster,errMessage,
        //	       "EnvelopeBunch time reset at z = %.3lf m with: %20.14le s, new offset: %20.14le s",
        //zAvg(),t,t+t_offset);
        writeError();
    }
    t_offset += new_dt;
    t        -= new_dt;
    return new_dt;
}

/** derivs for RK routine
 * Definition of equation set:
 *  ==========================
 Y[SLI_z]    = z      dz/dt  = beta*c*cos(a)
 cos(a) = 1 - (px0^2 + py0^2)/c2
 Y[SLI_beta] = beta   db/dt  = (e0/mc)*(E_acc + E_sc)/gamma^3
 Y[SLI_x]    = x      dx/dt  = px = Y[SLI_px]
 Y[SLI_px]   = px     dpx/dt = f(x,beta) - (beta*gamma^2(db/dt)*px + Kr*x) 
 Y[SLI_y]    = y      dy/dt  = py = Y[SLI_py]
 Y[SLI_py]   = py     dpy/dt = f(y,beta) - (beta*gamma^2(db/dt)*py + Kr*y)
 Y[SLI_x0]   = x0     dx0/dt = px0 = Y[SLI_px0]
 Y[SLI_px0]  = px0    dpx0/dt= -(beta*gamma^2(db/dt)*px0) + Kt*x0
 Y[SLI_y0]   = y0     dy0/dt = py0 = Y[SLI_py0]
 Y[SLI_py0]  = py0    dpy0/dt= -(beta*gamma^2(db/dt)*py0) + Kt*y0

 space charge blowup:
 f(x,beta) = c^2*I/(2Ia)/(x*beta*gamma^3)
 */

void EnvelopeBunch::derivs(double tc, double Y[], double dYdt[])
{
    double g2    = 1.0/(1.0-pow(Y[SLI_beta],2));
    double g     = sqrt(g2);
    double ecbmg = -1.0*Physics::q_e*Physics::c*Y[SLI_beta]/(g*Physics::EMASS);
    double c2    = Physics::c*Physics::c;

    /// minimum spot-size due to emittance
    double enxc2 = pow((emtnx0 + emtbx0)*Physics::c/(Y[SLI_beta]*g),2);
    double enyc2 = pow((emtny0 + emtby0)*Physics::c/(Y[SLI_beta]*g),2);

    /// transverse space charge
    /// somewhat strange: I expected: c*c*I/(2*Ia) (R. Bakker)
    ///  kpc  = 0.375*c*c*I->max()/Ia;
    double kpc  = 0.5*Physics::c*Physics::c*I->max()/Physics::Ia;

    // dYdt[SLI_z] = Y[SLI_beta]*c*sqrt(1.0 - (pow(Y[SLI_px0],2) + pow(Y[SLI_py0],2))/pow(c*Y[SLI_beta],2));
    //dYdt[SLI_z]    = Y[SLI_beta]*c*cos(Y[SLI_px0]/Y[SLI_beta]/c)*cos(Y[SLI_py0]/Y[SLI_beta]/c);
    dYdt[SLI_z]    = Y[SLI_beta]*Physics::c*cos(sqrt(pow(Y[SLI_px0],2)+pow(Y[SLI_py0],2))/Y[SLI_beta]/Physics::c);

    //update beta
    //FIXME: -Esl(2)!!!!
    //FIXME: do this here?? TW same probem??
    dYdt[SLI_beta] = Physics::e0mc*(-Esl(2) + Esct[cS] + Ezw[cS])/pow(g,3);

    /// beta * gamma^2 * dbeta/dt
    double bg2dbdt = Y[SLI_beta]*g2*dYdt[SLI_beta];

    if(solver & sv_radial) {
        dYdt[SLI_x]  = Y[SLI_px];
        dYdt[SLI_px] = (kpc*G[cS]/(Y[SLI_x]*Y[SLI_beta]*pow(g,3))) +
                       (enxc2/pow(Y[SLI_x],3)) - (KRsl(0)*Y[SLI_x]) - (bg2dbdt*Y[SLI_px]);
        dYdt[SLI_y]  = Y[SLI_py];
        dYdt[SLI_py] = (kpc*G[cS]/(Y[SLI_y]*Y[SLI_beta]*pow(g,3))) +
                       (enyc2/pow(Y[SLI_y],3)) - (KRsl(1)*Y[SLI_y]) - (bg2dbdt*Y[SLI_py]);
    } else {
        dYdt[SLI_x]  = Y[SLI_px];
        dYdt[SLI_px] = 0.0;
        dYdt[SLI_y]  = Y[SLI_py];
        dYdt[SLI_py] = 0.0;
    }

    if(solver & sv_offaxis) {
        dYdt[SLI_x0]  = Y[SLI_px0];
        dYdt[SLI_px0] = -KTsl(0) - (bg2dbdt*Y[SLI_px0]) + Physics::e0m*(g*Exw[cS]);
        dYdt[SLI_y0]  = Y[SLI_py0];
        dYdt[SLI_py0] = -KTsl(1) - (bg2dbdt*Y[SLI_py0]) + Physics::e0m*(g*Eyw[cS]);
    } else {
        dYdt[SLI_x0]  = Y[SLI_px0];
        dYdt[SLI_px0] = 0.0;
        dYdt[SLI_y0]  = Y[SLI_py0];
        dYdt[SLI_py0] = 0.0;
    }
}

/** Define a function to set parameters of each slice, where:
 *  bsX = width in X direction
 *  bsY = width in Y direction
 */
//ff
void EnvelopeBunch::DefineSlices(int no, double beta, double pos, double bsX, double bsY) 
{
    s[no]->p[1]=beta;
    s[no]->p[0]=pos;
    s[no]->p[2]=bsX;
    s[no]->p[4]=bsY;
}

void EnvelopeBunch::writeSlices() 
{
    int slices;
    double temp;
    Inform msg("writeSlices");

    slices=getN();
    msg << "Write output for " << slices << " slices:" << endl;
    msg << "\tn\tx\ty\tz-coords\tz-width" << endl;
    temp = s[0]->p[SLI_z];

    for (int i=0;i<slices;i++) {
        msg << "\t" << i << "\t" << s[i]->p[SLI_x] << "\t" << s[i]->p[SLI_y] << "\t" << s[i]->p[SLI_z] << "\t";
        if (i>0) {
            temp=s[i]->p[SLI_z]-temp;
            msg << temp << endl;
            temp=s[i]->p[SLI_z];
        } else {
            msg << "-" << endl;
        }
    }
    msg << "Finished writeSlice" << endl;
}

/** bunch of set and get functions for SLPartbunch 
*/

void EnvelopeBunch::setZ(int i, double coo) 
{
    s[i]->p[SLI_z] = coo;
}

double EnvelopeBunch::getZ(int i) 
{
    Inform msg("getZ(i)");
    double zCoord=0;
    if (i<getN()) {
        zCoord=s[i]->p[SLI_z];
        return zCoord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return zCoord;
    }
}

//ff
double EnvelopeBunch::getX(int i) 
{
    Inform msg("getX(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_x];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}

//ff
double EnvelopeBunch::getY(int i) 
{
    Inform msg("getY(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_y];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}

//ff
double EnvelopeBunch::getX0(int i) 
{
    Inform msg("getX(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_x0];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}

//ff
double EnvelopeBunch::getY0(int i) 
{
    Inform msg("getY(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_y0];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}


double EnvelopeBunch::getPx(int i) 
{
    Inform msg("getPx(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_px];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}

double EnvelopeBunch::getPy(int i) 
{
    Inform msg("getPy(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_py];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}

double EnvelopeBunch::getPz(int i) 
{
    Inform msg("getPz(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_beta]*Physics::m_e*s[i]->gamma();
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}

double EnvelopeBunch::getPx0(int i) 
{
    Inform msg("getPz(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_px0];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}

double EnvelopeBunch::getPy0(int i) 
{
    Inform msg("getPz(i)");
    double Coord=0;
    if (i<getN()) {
        Coord=s[i]->p[SLI_py0];
        return Coord;
    } else {
        msg << "ERROR: Not enough slices!" << endl;
        return Coord;
    }
}


void EnvelopeBunch::updateFields() 
{
    Inform msg("updateFields()");
    static unsigned long nCalled = 0;
    double Esc=0;

    // Calculate the current profile

    if (Q > 0.0) {

        calcI(); 

        /** the following assumes space-charges do not change significantly
         *  over nSc steps 
         */

        cSpaceCharge();

        /** accordingly also the wake-fiels are considered static over a 
         *  single step 
         */				
        //    cWake();

        if (n > 40) {
            /** smooth Esct to get rid of numerical noise
            */
            //sgSmooth(Esct,n,n/20,n/20,0,4);
        }
    }

    ++nCalled;
} 

/** alternative timestep-function DOES NOT WORK 
*/

/*
   void EnvelopeBunch::runSS(int i, double _zcat, double tStep) {
   static unsigned long 
   nCalled = 0;
   static int 
   msgParsed = 0;

   int
   nok,nbad,
   i,iStart,i0,di;
   double
   eps    = 1.0e-4,        // default accuricy integration 
   dt;                     // integration step-size (tStep by default)
   Slice 
 *sp;

// mark this EnvelopeBunch & the Elements globally
itsBunch = this;
// fE       = fe;
zCat     = _zCat;

backup(); // backup last stage before new execution

cS = i;    // make the current slice index global in this class
sp = s[i];

//Assign values of K for certain slices

KRsl=KR[i];
KTsl=KT[i];
Esl=EF[i];
Bsl=BF[i];

int    ode_result;
double epsLocal;

epsLocal   = eps;        // set default allowed error for integration
ode_result = 1;          // mark that the integration was not succesfull yet

while (ode_result == 1) {

if (solver & sv_fixedStep) {
rk4(sp->p,SLNPAR,t,dt,Gderivs);
ode_result = 0;
} else {
ode_result = odeint(sp->p,SLNPAR,
t,t+dt,epsLocal,0.1*dt,0.0,
&nok,&nbad,Gderivs);
}

if (ode_result != 0) {
sp->restore();     // restore the backup
epsLocal *= 10.0;
}
}
if (ode_result == 1) { // use fixed step integration if dynamic fails
rk4(sp->p,SLNPAR,t,dt,Gderivs);

if (msgParsed < 2) {
//writeError(errModeAll,errWarning,
//	     "EnvelopeBunch::run() Switched to fixed step RK rountine for solving of DE at slice %d at z = %lf m, nC=%d (g = %.10le). ONLY FIRST OCCURENCE MARKED!",
//	     i,sp->p[SLI_z],nCalled+1,sp->gamma());
writeError();
msgParsed = 2;
}
} else if ((epsLocal != eps) && (msgParsed == 0)) {
//writeError(errModeAll,errWarning,
//	     "EnvelopeBunch::run() integration accuracy relaxed to %le for slice %d at z = %lf m, nC=%d. ONLY FIRST OCCURENCE MARKED!",
//	     epsLocal,i,sp->p[SLI_z],nCalled+1);
writeError();
msgParsed = 1;
}

if (s[i]->check()) {
    //writeError(errModeMaster,errWarning,
    //"Slice %d no longer valid at z = %.3lf (beta = %.5lf)",
    //i,s[i]->p_old[SLI_z],s[i]->p_old[SLI_beta]);
    writeError(); 
};

// mark that slices might not be synchronized (and space charge accordingly)
dStat &= (!(ds_slicesSynchronized | ds_spaceCharge));

// mark calling of this function + update vars
t += tStep; // update time
++nCalled;

// subtract average orbit for when tracking along the s-axis
if (solver & sv_s_path) {
    int
        nV = 0;
    double 
        ga   = 0.0,
             x0a  = 0.0, px0a = 0.0,
             y0a  = 0.0, py0a = 0.0;
    double
        beta,fi_x,fi_y;

    // calc average over 80% of the bunch (ignore head-tail effects)
    for (i=((int) (0.1*n)); i < ((int) (0.9*n)); i++) {
        sp  = s[i];
        if (/*(sp->p[SLI_z] >= zCat) && sp->valid) {
              ++nV;
              ga  += sp->gamma();
              x0a += sp->p[SLI_x0]; px0a += sp->p[SLI_px0];
              y0a += sp->p[SLI_y0]; py0a += sp->p[SLI_py0];
              }
              }
              if (nV > 0) {
              ga  = ga/nV;
              x0a = x0a/nV; px0a = px0a/nV;
              y0a = y0a/nV; py0a = py0a/nV;
              } else {
            //writeError(errModeAll,errWarning,
            //		 "EnvelopeBunch::run() No valid slices to subtract average");
            writeError();
            }
            beta = sqrt(1.0 - (1.0/pow(ga,2)));
            fi_x = px0a/c/beta;
            fi_y = py0a/c/beta;

            dx0 -= x0a; dfi_x -= fi_x;
            dy0 -= y0a; dfi_y -= fi_y;
            for (i=0; i<n; i++) {
            sp = s[i];

            sp->p[SLI_x0]  -= x0a; sp->p[SLI_px0] -= px0a;
            sp->p[SLI_y0]  -= y0a; sp->p[SLI_py0] -= py0a;

            sp->p[SLI_z]   += (sp->p[SLI_x0]*sin(fi_x) + 
            sp->p[SLI_y0]*sin(fi_y));
            }
            }

            }*/


/** returns time of the bunch
 *  \return time of the bunch
*/
double EnvelopeBunch::getT() 
{
    return t;
}

/** 
 * \brief The main timestep function
 * gathers all parameters and integrates. Then writes new parameter in the bunch.
*/
void EnvelopeBunch::run(double _zCat, double tStep)
{
    Inform msg("run");
    static unsigned long nCalled = 0;
    static int msgParsed = 0;

    int nok,nbad,i,iStart,i0,di;
    // default accuracy of integration 
    double eps = 1.0e-4;
    // integration step-size (tStep by default)
    double dt;
    EnvelopeSlice *sp;

    itsBunch = this;
    // fE       = fe;
    zCat     = _zCat;

    i0 = 0;
    di = 1;

    backup(); // backup last stage before new execution

/*
    /// modify stepsize in the case slices are exiting the cathode

    dt     = tStep;
    iStart = -1;
    // check if we still have to do emission
    if (s[0]->p[SLI_z] < zCat) {
        msg << "slice 0 behind cathode" << endl;
        // find last slice still behind cathode
        i = n-1;
        while((i>0) && (s[i]->p[SLI_z] >= zCat)) 
            --i;

        // emit first slice
        if(i >= 0) {
            iStart = i;
            // find out if closest slice will be emitted in this time step
            dt = (zCat - s[i]->p[SLI_z])/(s[i]->p[SLI_beta]*Physics::c);
            if (dt > tStep) {
                dt = tStep;
                // normal advance first slice behind cathode with tStep if not
                // emitted in this step
                s[iStart]->p[SLI_z] += (Physics::c*dt*s[i]->p[SLI_beta]);
            } else {
                // if slice would be emitted in this time step, set z = zCat
                s[iStart]->p[SLI_z]  = zCat;
            }
            // step all other slices still behind cathode (will not be emitted
            // in this time step anyway)
            for (i=0; i<iStart; i++) {
                s[i]->p[SLI_z] += (Physics::c*dt*s[i]->p[SLI_beta]);
            }
        }
    }
*/    
    
    dt = tStep;
    iStart = -1;
    // check if we (still) have to do emission
    if(s[0]->p[SLI_z] < zCat) {
        // find last slice still behind cathode
        i = n-1;
        while((i>0) && (s[i]->p[SLI_z] >= zCat)) 
            --i;

        // emit next slice
        if(i >= 0) {
            iStart = i;

            // find out if closest slice will be emitted in this time step
            dt = (zCat - s[i]->p[SLI_z])/(s[i]->p[SLI_beta]*Physics::c);

            // move next slide to cathode
            s[i]->p[SLI_z] = zCat;

            // step all other slices still behind cathode (will not be emitted
            // in this time step anyway)
            for(i=0; i<iStart; i++) {
                s[i]->p[SLI_z] += (Physics::c*dt*s[i]->p[SLI_beta]);
            }
        }
    } 

    for (i=0; i<n; i++) {
        // make the current slice index global in this class
        cS = i;   
        sp = s[i];

        // Assign values of K for certain slices
        KRsl = KR[i];
        KTsl = KT[i];
        Esl = EF[i];
        Bsl = BF[i];

        // only for slices already emitted
        if(i > iStart) {
            // set default allowed error for integration
            double epsLocal = eps;
            // mark that the integration was not successful yet
            int ode_result = 1;

            while(ode_result == 1) {

                if(solver & sv_fixedStep) {
                    rk4(sp->p,SLNPAR,t,dt,Gderivs);
                    ode_result = 0;
                } else {
                    ode_result = odeint(sp->p,SLNPAR, t,t+dt,epsLocal,0.1*dt,0.0, &nok,&nbad,Gderivs);
                }

                if(ode_result != 0) {
                    // restore the backup
                    sp->restore();     
                    epsLocal *= 10.0;
                }
            }

            if(ode_result == 1) {
                // use fixed step integration if dynamic fails
                rk4(sp->p,SLNPAR,t,dt,Gderivs);

                if (msgParsed < 2) {
                    msg << "EnvelopeBunch::run() Switched to fixed step RK rountine for solving of DE at slice " << i << endl;
                    msgParsed = 2;
                }
            } else if((epsLocal != eps) && (msgParsed == 0)) {
                writeError();
                msg << "EnvelopeBunch::run() integration accuracy relaxed to " << epsLocal << " for slice " << i << "ONLY FIRST OCCURENCE MARKED!" << endl;
                msgParsed = 1;
            }
        }

        if (s[i]->check()) {
            msg << "Slice " << i << " no longer valid at z = " <<  s[i]->p_old[SLI_z] << endl;
        }
    }
    // mark that slices might not be synchronized (and space charge accordingly)
    dStat &= (!(ds_slicesSynchronized | ds_spaceCharge));

    /// mark calling of this function + update vars
    t += tStep; // update time
    ++nCalled;

    /// subtract average orbit for when tracking along the s-axis
    if (solver & sv_s_path) {
        int nV = 0;
        double ga = 0.0, x0a = 0.0, px0a = 0.0, y0a = 0.0, py0a = 0.0;
        double beta,fi_x,fi_y;

        //FIXME: BET calculates only 80 %, OPAL doesn't ?

        for (i=0;i<n;i++) {
            sp  = s[i];
            if ((sp->p[SLI_z] >= zCat) && sp->valid) {
                ++nV;
                ga  += sp->gamma();
                x0a += sp->p[SLI_x0];
                y0a += sp->p[SLI_y0];
                px0a += sp->p[SLI_px0];
                py0a += sp->p[SLI_py0];
            }
        }
        if (nV > 0) {
            ga  = ga/nV;
            x0a = x0a/nV;
            px0a = px0a/nV;
            y0a = y0a/nV;
            py0a = py0a/nV;
        } else {
            writeError();
            msg << "EnvelopeBunch::run() No valid slices to subtract average" << endl;
        }
        beta = sqrt(1.0 - (1.0/pow(ga,2)));
        fi_x = px0a/Physics::c/beta;
        fi_y = py0a/Physics::c/beta;

        dx0 -= x0a;
        dy0 -= y0a;
        dfi_x -= fi_x;
        dfi_y -= fi_y;
        for (i=0; i<n; i++) {
            sp = s[i];

            sp->p[SLI_x0] -= x0a;
            sp->p[SLI_y0] -= y0a;
            sp->p[SLI_px0] -= px0a;
            sp->p[SLI_py0] -= py0a;
            sp->p[SLI_z] += (sp->p[SLI_x0]*sin(fi_x) + sp->p[SLI_y0]*sin(fi_y));
        }
    }
}

