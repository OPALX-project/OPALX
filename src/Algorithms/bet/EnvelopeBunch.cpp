#include <iostream>
#include <cfloat>
#include <fstream>
#include <cmath>
#include <string>
#include <assert.h>
//FIXME: replace with IPPL Vector_t
#include <vector>
//FIXME: remove
#include <mpi.h>

#include "Algorithms/bet/math/root.h"     // root finding routines
#include "Algorithms/bet/math/sort.h"     // sorting routines
#include "Algorithms/bet/math/linfit.h"   // linear fitting routines
#include "Algorithms/bet/math/savgol.h"   // savgol smoothing routine
#include "Algorithms/bet/math/rk.h"       // Runge-Kutta Integration

#include "Algorithms/bet/EnvelopeBunch.h"
#include <Physics/Physics.h>

extern Inform *gmsg;

/** for space charge
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
    beta = 0.85 -> 460 keV
*/
#define BETA_MIN1 0.30     // minimum beta-value for space-charge calculations: start
#define BETA_MIN2 0.45     // minimum beta-value for space-charge calculations: full impact

// Hack allows odeint in rk.C to be called with a class member function
static EnvelopeBunch *thisBunch = NULL;  // pointer to access calling bunch
static void Gderivs(double t, double Y[], double dYdt[]) { thisBunch->derivs(t, Y, dYdt); }
static double Gcur(double z) { return thisBunch->I->get(z, itype_lin); }
static double rootValue = 0.0;

// used in setLShape for Gaussian
static void erfRoot(double x, double *fn, double *df) {
    double v = erfc(fabs(x));
    double eps = 1.0e-05;

    *fn = v - rootValue;
    *df = (erfc(fabs(x) + eps) - v) / eps;
}

void EnvelopeBunch::calcBeamParameters() {
    Inform msg("calcBeamParameters");
    IpplTimings::startTimer(statParamTimer_m);
    double ex, ey, nex, ney, b0, bRms, bMax, bMin, g0, dgdt, gInc;
    double RxRms, RyRms, Px, PxRms, PxMax, PxMin, Py, PyRms, PyMax, PyMin;
    double Pz, PzMax, PzMin, PzRms;
    double x0Rms, y0Rms, zRms, zMax, zMin, I0, IRms, IMax, IMin;

    runStats(sp_beta, &b0, &bMax, &bMin, &bRms, &nValid_m);
    runStats(sp_I, &I0, &IMax, &IMin, &IRms, &nValid_m);
    runStats(sp_z, &z0_m, &zMax, &zMin, &zRms, &nValid_m);
    runStats(sp_Pz, &Pz, &PzMax, &PzMin, &PzRms, &nValid_m);

    if(solver & sv_radial) {
        runStats(sp_Rx, &Rx_m, &RxMax_m, &RxMin_m, &RxRms, &nValid_m);
        runStats(sp_Ry, &Ry_m, &RyMax_m, &RyMin_m, &RyRms, &nValid_m);
        runStats(sp_Px, &Px, &PxMax, &PxMin, &PxRms, &nValid_m);
        runStats(sp_Py, &Py, &PyMax, &PyMin, &PyRms, &nValid_m);
        calcEmittance(&nex, &ney, &ex, &ey, &nValid_m);
    }

    if(solver & sv_offaxis) {
        runStats(sp_x0, &x0_m, &x0Max_m, &x0Min_m, &x0Rms, &nValid_m);
        runStats(sp_y0, &y0_m, &y0Max_m, &y0Min_m, &y0Rms, &nValid_m);
    }

    calcEnergyChirp(&g0, &dgdt, &gInc, &nValid_m);
    double Bfz = AvBField();
    double Efz = AvEField();
    double mc2e = 1.0e-6 * Physics::EMASS * Physics::c * Physics::c / Physics::q_e;

    E_m = mc2e * (g0 - 1.0);
    dEdt_m = 1.0e-12 * mc2e * dgdt;
    Einc_m = mc2e * gInc;
    tau_m =  zRms / Physics::c;
    I_m = IMax;
    Irms_m = Q_m * nValid_m * Physics::c / (zRms * sqrt(Physics::two_pi) * numSlices_m);
    Px_m = Px / Physics::c;
    Py_m = Py / Physics::c;
    dx0_m = dx0;
    dy0_m = dy0;
    dfi_x_m = 180.0 * dfi_x / Physics::pi;
    dfi_y_m = 180.0 * dfi_y / Physics::pi;
    Ez_m = 1.0e-6 * Efz;
    Bz_m = Bfz;

    //in [mrad]
    emtn_m = Vector_t(ex, ey, 0.0);
    norm_emtn_m = Vector_t(nex, ney, 0.0);

    maxX_m = Vector_t(RxMax_m, RyMax_m, zMax);
    maxP_m = Vector_t(PxMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, PyMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, PzMax);

    //minX in the T-T sense is -RxMax_m
    minX_m = Vector_t(-RxMax_m, -RyMax_m, zMin);
    minP_m = Vector_t(-PxMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, -PyMax * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi, PzMin);
    //minP_m = Vector_t(-maxP_m[0], -maxP_m[1], PzMin);

    /* PxRms is the rms of the divergence in x direction in [rad]: x'
     * x' = Px/P0 -> Px = x' * P0 = x' * \beta/c*E (E is the particle total
     * energy)
     * Pz* in [\beta\gamma]
     * \pi from [rad]
     */
    sigmax_m = Vector_t(RxRms / 2.0, RyRms / 2.0, zRms);
    sigmap_m = Vector_t(PxRms * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi / 2.0, PyRms * E_m * sqrt((g0 + 1.0) / (g0 - 1.0)) / Physics::c * Physics::pi / 2.0, PzRms);

    IpplTimings::stopTimer(statParamTimer_m);
}

void EnvelopeBunch::runStats(EnvelopeBunchParameter sp, double *xAvg, double *xMax, double *xMin, double *rms, int *nValid) {
    int i, nV = 0;
    double *v = vector(numMySlices_m);

    //FIXME: why from 1 to n-1??
    switch(sp) {
        case sp_beta:      // normalized velocity (total) [-]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_beta];
            }
            break;
        case sp_gamma:     // Lorenz factor
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->gamma();
            }
            break;
        case sp_z:         // slice position [m]
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_z];
            }
            break;
        case sp_I:         // slice position [m]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_beta] > BETA_MIN1) && ((s[i]->p[SLI_z] > zCat) && s[i]->valid))
                    v[nV++] = I->get(s[i]->p[SLI_z], itype_lin);
            }
            break;
        case sp_Rx:        // beam size x [m]
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = 2.0 * s[i]->p[SLI_x];
            }
            break;
        case sp_Ry:        // beam size y [m]
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = 2.0 * s[i]->p[SLI_y];
            }
            break;
        case sp_Px:        // beam divergence x
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_px];
            }
            break;
        case sp_Py:        // beam divergence y
            for(i = 0; i < numMySlices_m - 0; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_py];
            }
            break;
        case sp_Pz:
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_beta] * s[i]->gamma();
            }
            break;
        case sp_x0:        // position centroid x [m]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_x0];
            }
            break;
        case sp_y0:        // position centroid y [m]
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_y0];
            }
            break;
        case sp_px0:       // angular deflection centroid x
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_px0];
            }
            break;
        case sp_py0:      // angular deflection centroid y
            for(i = 1; i < numMySlices_m - 1; i++) {
                if((s[i]->p[SLI_z] > zCat) && s[i]->valid) v[nV++] = s[i]->p[SLI_py0];
            }
            break;
        default :
            throw OpalException("EnvelopeBunch", "EnvelopeBunch::runStats() Undefined label (programming error)");
            break;
    }

    int nVTot = nV;
    reduce(nVTot, nVTot, OpAddAssign());
    if(nVTot <= 0) {
        *xAvg = 0.0;
        *xMax = 0.0;
        *xMin = 0.0;
        *rms = 0.0;
        *nValid = 0;
    } else {
        double M1 = 0.0;
        double M2 = 0.0;
        double maxv = numeric_limits<double>::min();
        double minv = numeric_limits<double>::max();
        if(nV > 0) {
            M1 = v[0];
            M2 = v[0] * v[0];
            maxv = v[0];
            minv = v[0];
            for(i = 1; i < nV; i++) {
                M1 += v[i];
                M2 += v[i] * v[i];
                maxv = max(maxv, v[i]);
                minv = min(minv, v[i]);
            }
        }

        reduce(M1, M1, OpAddAssign());
        reduce(M2, M2, OpAddAssign());
        reduce(maxv, maxv, OpMaxAssign());
        reduce(minv, minv, OpMinAssign());

        *xAvg = M1 / nVTot;
        *xMax = maxv;
        *xMin = minv;

        //XXX: ok so this causes problems. E.g in the case of all transversal
        //components we cannot compare a rms_rx with rms_x (particles)
        //produced by the T-Tracker

        //in case of transversal stats we want to calculate the rms
        //*rms = sqrt(M2/nVTot);

        //else a sigma
        //*sigma = sqrt(M2/nVTot - M1*M1/(nVTot*nVTot));
        //this is sigma:
        //*rms = sqrt(M2/nVTot - M1*M1/(nVTot*nVTot));
        *nValid = nVTot;

        if(sp == sp_Rx ||
           sp == sp_Ry ||
           sp == sp_Px ||
           sp == sp_Py)
            *rms = sqrt(M2 / nVTot); //2.0: half of the particles
        else
            *rms = sqrt(M2 / nVTot - M1 * M1 / (nVTot * nVTot));
    }
    free(v);
}

void EnvelopeBunch::calcEmittance(double *emtnx, double *emtny, double *emtx, double *emty, int *nValid) {
    double sx = 0.0;
    double sxp = 0.0;
    double sxxp = 0.0;
    double sy = 0.0;
    double syp = 0.0;
    double syyp = 0.0;
    double betagamma = 0.0;

    // find the amount of active slices
    int nV = 0;
    for(int i = 0; i < numMySlices_m; i++) {
        if((s[i]->p[SLI_z] > zCat) && s[i]->valid)
            nV++;
    }

    if(nV > 0) {
        int i1 = nV;
        nV = 0;
        double bg = 0.0;
        for(int i = 0; i < i1; i++) {
            if((s[i]->p[SLI_z] > zCat) && s[i]->valid) {
                nV++;

                if(solver & sv_radial) {
                    assert(i < numMySlices_m);
                    EnvelopeSlice *sp = s[i];
                    bg = sp->p[SLI_beta] * sp->gamma();

                    double pbc = bg * sp->p[SLI_px] / (sp->p[SLI_beta] * Physics::c);
                    sx   += sp->p[SLI_x] * sp->p[SLI_x];
                    sxp  += pbc * pbc;
                    sxxp += sp->p[SLI_x] * pbc;

                    pbc   = bg * sp->p[SLI_py] / (sp->p[SLI_beta] * Physics::c);
                    sy   += sp->p[SLI_y] * sp->p[SLI_y];
                    syp  += pbc * pbc;
                    syyp += sp->p[SLI_y] * pbc;

                    betagamma += sqrt(1 + sp->p[SLI_px] * sp->p[SLI_px] + sp->p[SLI_py] * sp->p[SLI_py]);
                }
            }
        }
    }

    int nVToT = nV;
    reduce(nVToT, nVToT, OpAddAssign());
    if(nVToT == 0) {
        *emtnx = 0.0;
        *emtny = 0.0;
        *emtx = 0.0;
        *emty = 0.0;
        *nValid = 0;
    } else {
        reduce(sx, sx, OpAddAssign());
        reduce(sy, sy, OpAddAssign());
        reduce(sxp, sxp, OpAddAssign());
        reduce(syp, syp, OpAddAssign());
        reduce(sxxp, sxxp, OpAddAssign());
        reduce(syyp, syyp, OpAddAssign());
        reduce(betagamma, betagamma, OpAddAssign());
        sx /= nVToT;
        sy /= nVToT;
        sxp /= nVToT;
        syp /= nVToT;
        sxxp /= nVToT;
        syyp /= nVToT;

        *emtnx = sqrt(sx * sxp - sxxp * sxxp + emtnx0 * emtnx0 + emtbx0 * emtbx0);
        *emtny = sqrt(sy * syp - syyp * syyp + emtny0 * emtny0 + emtby0 * emtby0);

        betagamma /= nVToT;
        betagamma *= sqrt(1.0 - (1 / betagamma) * (1 / betagamma));
        *emtx = *emtnx / betagamma;
        *emty = *emtny / betagamma;
        *nValid = nVToT;
    }
}

void EnvelopeBunch::calcEnergyChirp(double *g0, double *dgdt, double *gInc, int *nValid) {
    double *dt = vector(numMySlices_m);  /// z --> time
    double *b  = vector(numMySlices_m);  /// beta
    double *g  = vector(numMySlices_m);  /// gamma

    // defaults
    *g0   = 1.0;
    *dgdt = 0.0;
    *gInc = 0.0;

    double zAvg = 0.0;
    double gAvg = 0.0;
    int j = 0;
    for(int i = 0; i < numMySlices_m; i++) {
        EnvelopeSlice *cs = s[i];
        if((cs->valid) && (cs->p[SLI_z] > zCat)) {
            zAvg    += cs->p[SLI_z];
            dt[i-j]  = cs->p[SLI_z];
            b[i-j]   = cs->p[SLI_beta];
            g[i-j]   = cs->gamma();
            gAvg    += g[i-j];
        } else
            ++j;
    }
    int nV = numMySlices_m - j;

    int nVTot = nV;
    reduce(nVTot, nVTot, OpAddAssign());
    if(nVTot > 0) {
        reduce(gAvg, gAvg, OpAddAssign());
        reduce(zAvg, zAvg, OpAddAssign());
        gAvg = gAvg / nVTot;
        zAvg = zAvg / nVTot;
        *g0  = gAvg;
    }

    // FIXME: working with global arrays
    // make dt, b and g global
    if(nVTot > 2) {
        double *dtG = vector(nVTot);  // z --> time
        double *bG  = vector(nVTot);  // beta
        double *gG  = vector(nVTot);  // gamma

        int numproc = Ippl::Comm->getNodes();
        int *numsend = (int *)malloc(numproc * sizeof(int));
        int *offsets = (int *)malloc(numproc * sizeof(int));
        int *offsetsG = (int *)malloc(numproc * sizeof(int));
        int *zeros = (int *)malloc(numproc * sizeof(int));

        for(int i = 0; i < numproc; i++) {
            zeros[i] = 0;
            numsend[i] = nV;
        }

        MPI_Allgather(&nV, 1, MPI_INT, &offsets[0], 1, MPI_INT, MPI_COMM_WORLD);
        offsetsG[0] = 0;
        for(int i = 1; i < numproc; i++) {
            if(offsets[i-1] == 0)
                offsetsG[i] = 0;
            else
                offsetsG[i] = offsetsG[i-1] + offsets[i-1];
        }

        MPI_Alltoallv(&dt[0], numsend, zeros, MPI_DOUBLE, &dtG[0], offsets, offsetsG, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(&b[0], numsend, zeros, MPI_DOUBLE, &bG[0], offsets, offsetsG, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(&g[0], numsend, zeros, MPI_DOUBLE, &gG[0], offsets, offsetsG, MPI_DOUBLE, MPI_COMM_WORLD);

        free(offsets);
        free(offsetsG);
        free(numsend);
        free(zeros);

        double dum1, dum2, dum3, dum4, rms, gZero, gt;

        // convert z to t
        for(int i = 0; i < nVTot; i++) {
            dtG[i] = (dtG[i] - zAvg) / (bG[i] * Physics::c);
        }

        // chrip and uncorrelated energy sread
        linfit(&dtG[0], &gG[0], nVTot, &gZero, &gt, &dum2, &dum3, &dum4);
        *dgdt = gt;

        rms = 0.0;
        for(int i = 0; i < nVTot; i++) {
            rms += pow(gG[i] - gZero - gt * dtG[i], 2);
        }

        *gInc = sqrt(rms / nVTot);

        free(dtG);
        free(bG);
        free(gG);
    }

    free(g);
    free(b);
    free(dt);
}

EnvelopeBunch::EnvelopeBunch(const PartData *ref):
    reference(ref),
    PartBunch(ref)
{}


EnvelopeBunch::EnvelopeBunch(const EnvelopeBunch &rhs):
    reference(rhs.reference),
    PartBunch(rhs.reference)
{}


EnvelopeBunch::EnvelopeBunch(const std::vector<Particle> &rhs, const PartData *ref):
    reference(ref),
    PartBunch(ref)
{}


EnvelopeBunch::~EnvelopeBunch() {
    delete[] LastSection;
    delete[] dt;
    for(int i = 0; i < numMySlices_m; i++)
        delete s[i];
    delete s;
    numMySlices_m = 0;
    numSlices_m = 0;

    if(Exw)  free(Exw);
    if(Eyw)  free(Eyw);
    if(Ezw)  free(Ezw);
    if(Esct) free(Esct);
    if(G)    free(G);
    if(I)    delete I;
    if(z_m)  free(z_m);
    if(b_m)  free(b_m);
}

void EnvelopeBunch::createSlices(int nSlice) {
    numSlices_m = nSlice;
    int rank = Ippl::Comm->myNode();
    int numproc = Ippl::Comm->getNodes();
    numMySlices_m = nSlice / numproc;
    if(numMySlices_m < 13) {
        if(rank == 0) {
            numMySlices_m = 14;
        } else {
            numMySlices_m = (nSlice - 14) / (numproc - 1);
            if(rank - 1 < (nSlice - 14) % (numproc - 1))
                numMySlices_m++;
        }
    } else {
        if(rank < nSlice % numproc)
            numMySlices_m++;
    }

    mySliceStartOffset_m = rank * ((int)numSlices_m / numproc);
    if(rank < numSlices_m % numproc)
        mySliceStartOffset_m += rank;
    else
        mySliceStartOffset_m += numSlices_m % numproc;
    mySliceEndOffset_m = mySliceStartOffset_m + numMySlices_m - 1;

    LastSection = new long [getLocalNum()];
    dt = new double [getLocalNum()];

    KR = new Vector_t[numMySlices_m];
    KT = new Vector_t[numMySlices_m];
    EF = new Vector_t[numMySlices_m];
    BF = new Vector_t[numMySlices_m];

    z_m = vector(numSlices_m);
    b_m = vector(numSlices_m);

    for(int i = 0; i < numMySlices_m; i++) {
        KR[i] = (0, 0, 0);
        KT[i] = (0, 0, 0);
        BF[i] = (0, 0, 0);
        EF[i] = (0, 0, 0);
    }

    // set default DE solver method
    solver = sv_radial | sv_offaxis | sv_lwakes | sv_twakes;

    if(numSlices_m < 14)
        throw OpalException("EnvelopeBunch::createSlices", "use more than 13 slices");

    // defaults:
    Q_m        = 0.0;                   // no charge
    t        = 0.0;                   // t = 0 s
    t_offset = 0.0;                   // offset time by tReset function
    emtnx0   = 0.0;
    emtny0 = 0.0;     // no intrinsic emittance
    emtbx0   = 0.0;
    emtby0 = 0.0;     // no intrinsic emittance Bush effect
    Bz0      = 0.0;                   // no magnetic field on creation (cathode)

    dx0      = 0.0;
    dy0    = 0.0;     // no offset of coordinate system
    dfi_x    = 0.0;
    dfi_y  = 0.0;     // no rotation of coordinate system

    s = new EnvelopeSliceP[numMySlices_m];
    for(int i = 0; i < numMySlices_m; i++)
        s[i] = new EnvelopeSlice();

    /// allocate memory for internal arrays
    Esct = vector(numMySlices_m);
    G    = vector(numMySlices_m);
    Exw  = vector(numMySlices_m);
    Eyw  = vector(numMySlices_m);
    Ezw  = vector(numMySlices_m);
    for(int i = 0; i < numMySlices_m; i++) {
        Esct[i] = 0.0;
        G[i]   = 0.0;
        Ezw[i]  = 0.0;
        Exw[i]  = 0.0;
        Eyw[i] = 0.0;
    }

    I     = NULL;
    I0avg = 0.0;
    dStat = ds_fieldsSynchronized | ds_slicesSynchronized;
}

void EnvelopeBunch::setLShape(EnvelopeBunchShape shape, double z0, double w, double frac) {
    Inform msg("setLShape");
    int n2 = numSlices_m / 2;
    double sqr2 = sqrt(2.0), dz0, sl, v, vMax;
    double ffact = 1.0;
    double dztot;
    switch(shape) {
        case bsRect:

#ifndef BETEMISSIONMODEL
            dztot = Physics::c * w * s[0]->p[SLI_beta];
            for(int i = 0; i < numMySlices_m; i++) {
                s[i]->p[SLI_z] = ((numSlices_m - (mySliceStartOffset_m + i)) * dztot) / numSlices_m;
                //s[i]->p[SLI_z] = -1.0;
            }
            I0avg = Q_m * (w < 0.0 ? 1.0 : Physics::c) / fabs(2.0 * w);
            break;
#else
            //XXX: old bet emission model only serial at the moment
            sl   = fabs(frac) > fabs(w) ? w / 2.0 : frac;
            ffact = w < 0.0 ? Physics::c : 1.0;
            vMax = ffact * (fabs(w) - fabs(sl)) / 2.0;
            dz0  = ffact * (fabs(w) + 0.95 * fabs(sl)) / (numSlices_m - 1);
            sl   = ffact * fabs(sl);
            v    = 0.0;
            s[n2]->p[SLI_z] = z0;

            for(int i = 1; i <= numSlices_m / 2.0; i++) {
                if(v <= vMax)
                    v += dz0;
                else
                    v += (2.0 * dz0 / (1.0 + cos(Physics::u_two_pi * (v - vMax) / sl)));

                if(n2 + i <  numSlices_m)
                    s[n2+i]->p[SLI_z] = z0 + v * s[n2+i]->p[SLI_beta];
                if(n2 - i >= 0)
                    s[n2-i]->p[SLI_z] = z0 - v * s[n2-i]->p[SLI_beta];
            }

            /*
            sl   = fabs(frac) > fabs(w) ? w/2.0 : frac;
            ffact = w < 0.0 ? Physics::c : 1.0;
            vMax = ffact*(fabs(w) - fabs(sl))/2.0;
            dz0  = ffact*(fabs(w) + 0.95*fabs(sl))/(numSlices_m-1);
            sl   = ffact*fabs(sl);
            v    = 0.0;
            if(n2 >= mySliceStartOffset_m && n2 <= mySliceEndOffset_m) {
                //FIXME: why did Rene use z0 instead of z0*beta?
                //z0 *= s[n2-mySliceStartOffset_m]->p[SLI_beta];
                s[n2-mySliceStartOffset_m]->p[SLI_z] = z0;
            }

            for(int i=1; i<=numSlices_m/2; i++) {
                if(v <= vMax)
                    v += dz0;
                else
                    v += (2.0*dz0/(1.0+cos(Physics::u_two_pi*(v-vMax)/sl)));

                if((n2-i) >= mySliceStartOffset_m && (n2-i) <= mySliceEndOffset_m)
                    s[n2-i-mySliceStartOffset_m]->p[SLI_z] = z0 - v*s[n2-i-mySliceStartOffset_m]->p[SLI_beta];

                if((n2+i) >= mySliceStartOffset_m && (n2+i) <= mySliceEndOffset_m)
                    s[n2+i-mySliceStartOffset_m]->p[SLI_z] = z0 + v*s[n2+i-mySliceStartOffset_m]->p[SLI_beta];
            }
            */

            I0avg = Q_m * (w < 0.0 ? 1.0 : Physics::c) / fabs(2.0 * w);
            break;
#endif

        case bsGauss:
            if(n2 >= mySliceStartOffset_m && n2 <= mySliceEndOffset_m)
                s[n2]->p[SLI_z] = z0;

            for(int i = 1; i <= numSlices_m / 2; i++) {
                rootValue = 1.0 - 2.0 * i * frac / (numMySlices_m + 1);
                v = fabs(w) * sqr2 * findRoot(erfRoot, 0.0, 5.0, 1.0e-5) * (w < 0.0 ? Physics::c : 1.0);

                if((n2 + i) >= mySliceStartOffset_m && (n2 + i) <= mySliceEndOffset_m)
                    s[n2+i-mySliceStartOffset_m]->p[SLI_z] = z0 + v * s[n2+i-mySliceStartOffset_m]->p[SLI_beta];

                if((n2 - i) >= mySliceStartOffset_m && (n2 - i) <= mySliceEndOffset_m)
                    s[n2-i-mySliceStartOffset_m]->p[SLI_z] = z0 - v * s[n2-i-mySliceStartOffset_m]->p[SLI_beta];
            }

            I0avg = 0.0;
            break;
    }

    // save new state
    backup();
}

void EnvelopeBunch::setTShape(double enx, double eny, double rx, double ry, double b0) {
    Inform msg("setTshape");
    msg << "set SLI_x to " << rx / 2.0 << endl;
    msg << "set SLI_y to " << ry / 2.0 << endl;

    // set emittances
    emtnx0 = enx;
    emtny0 = eny;
    emtbx0 = Physics::q_e * rx * rx * Bz0 / (8.0 * Physics::EMASS * Physics::c);
    emtby0 = Physics::q_e * ry * ry * Bz0 / (8.0 * Physics::EMASS * Physics::c);
    msg << "set emtbx0 to " << emtbx0 << endl;
    msg << "set emtby0 to " << emtby0 << endl;

    // SLI_x = radius
    for(int i = 0; i < numMySlices_m; i++) {
        s[i]->p[SLI_x] = rx / 2.0;
        s[i]->p[SLI_px] = 0.0;
        s[i]->p[SLI_y] = ry / 2.0;
        s[i]->p[SLI_py] = 0.0;
    }

    // save new state
    backup();
}

void EnvelopeBunch::setTOffset(double x0, double px0, double y0, double py0) {
    for(int i = 0; i < numMySlices_m; i++) {
        s[i]->p[SLI_x0] = x0;
        s[i]->p[SLI_px0] = px0;
        s[i]->p[SLI_y0] = y0;
        s[i]->p[SLI_py0] = py0;
    }
}

void EnvelopeBunch::setEnergy(double E0, double dE) {
    double g0    = 1.0 + (Physics::q_e * E0 / (Physics::EMASS * Physics::c * Physics::c));
    double dg    = fabs(dE) * Physics::q_e / (Physics::EMASS * Physics::c * Physics::c);
    double z0    = zAvg();

    for(int i = 0; i < numMySlices_m; i++) {
        double g = g0 + (s[i]->p[SLI_z] - z0) * dg;
        s[i]->p[SLI_beta] = sqrt(1.0 - (1.0 / (g * g)));
    }

    // save new states
    backup();
}

void EnvelopeBunch::setSolver(int s) {
    solver = s;
}

void EnvelopeBunch::synchronizeSlices() {

    for(int i = 0; i < numSlices_m; i++) {
        z_m[i] = 0.0;
        b_m[i] = 0.0;
    }
    for(int i = 0; i < numMySlices_m; i++) {
        b_m[mySliceStartOffset_m+i] = s[i]->p[SLI_beta];
        z_m[mySliceStartOffset_m+i] = s[i]->p[SLI_z];
    }

    reduce(&(z_m[0]), &(z_m[0]) + numSlices_m, &(z_m[0]), OpAddAssign());
    reduce(&(b_m[0]), &(b_m[0]) + numSlices_m, &(b_m[0]), OpAddAssign());
}

void EnvelopeBunch::calcI() {
    Inform msg("calcI ");

    static int already_called = 0;
    if((dStat & ds_currentCalculated) || (already_called && (Q_m <= 0.0)))
        return;
    already_called = 1;

    // delete the old profile
    if(I)
        delete I;

    double *z1 = vector(numSlices_m);
    double *b = vector(numSlices_m);
    double bSum = 0.0;
    double dz2Sum = 0.0;
    int n1 = 0;

    for(int i = 0; i < numSlices_m; i++) {
        z1[i] = 0.0;
        b[i] = 0.0;
    }
    for(int i = 0; i < numSlices_m; i++) {
        if(b_m[i] > 0.0) {
            b[n1] = b_m[i];
            z1[n1] = z_m[i];
            if(n1 > 0)
                dz2Sum += pow(z1[n1] - z1[n1-1], 2);
            bSum += b_m[i];
            n1++;
        }
    }

    if(n1 < 2) {
        free(z1);
        free(b);
        msg << "n1 (= " << n1 << ") < 2" << endl;
        throw OpalException("EnvelopeBunch", "Insufficient points to calculate the current (n1)");
        I = new Profile(0.0);
        return;
    }

    double sigma_dz = sqrt(dz2Sum / (n1 - 1));
    double beta = bSum / n1;
    //sort z1 according to beta's
    sort2(z1, b, n1);

    double q = Q_m > 0.0 ? Q_m / numSlices_m : Physics::q_e;
    double dz = 0.0;

    // 1. determine current from the slice distance
    double *I1 = vector(n1);
    I1[0] = I1[n1-1] = 0.0;

    //limit the max current to 5x the sigma value
    // to reduce noise problems
    double dzMin = 0.2 * sigma_dz;
    for(int i = 1; i < n1 - 1; i++) {
        int j = 1;
        do {
            dz = fabs(z1[i+j] - z1[i-j]);
            j++;
        } while((dz < dzMin * (j - 1)) && ((i + j) < n1) && ((i - j) >= 0));

        if((dz >= dzMin * (j - 1)) && ((i + j) < n1) && ((i - j) >= 0))
            I1[i] = 0.25 * q * Physics::c * (b[i+j] + b[i-j]) / (dz * (j - 1));
        else
            I1[i] = I1[i-1];
    }

    I1[0]    = I1[1];
    I1[n1-1] = I1[n1-2];

    //2. Remove points with identical z-value
    // and then smooth the current profile
    double zMin = zTail();
    double zMax = zHead();
    dz = (zMax - zMin) / numSlices_m; // create a window of the average slice distance
    double *z2 = vector(n1);
    double *I2 = vector(n1);
    double Mz1 = 0.0;
    double MI1 = 0.0;
    int np = 0;
    int j = 0;

    // first value
    while((j < n1) && ((z1[j] - z1[0]) <= dz)) {
        Mz1 += z1[j];
        MI1 += I1[j];
        ++j;
        ++np;
    }
    z2[0] = Mz1 / np;
    I2[0] = MI1 / np;

    // following values
    int k = 0;
    for(int i = 1; i < n1; i++) {
        // add new points
        int j = 0;
        while(((i + j) < n1) && ((z1[i+j] - z1[i]) <= dz)) {
            if((z1[i+j] - z1[i-1]) > dz) {
                Mz1 += z1[i+j];
                MI1 += I1[i+j];
                ++np;
            }
            ++j;
        }

        // remove obsolete points
        j = 1;
        while(((i - j) >= 0) && ((z1[i-1] - z1[i-j]) < dz)) {
            if((z1[i] - z1[i-j]) > dz) {
                Mz1 -= z1[i-j];
                MI1 -= I1[i-j];
                --np;
            }
            ++j;
        }
        z2[i-k] = Mz1 / np;
        I2[i-k] = MI1 / np;

        // make sure there are no duplicate z coordinates
        if(z2[i-k] <= z2[i-k-1]) {
            I2[i-k-1] = 0.5 * (I2[i-k] + I2[i-k-1]);
            ++k;
        }
    }

    int n2 = n1 - k;
    if(n2 < 1) {
        *gmsg << "Insufficient points to calculate the current (m = " << n2 << ")" << endl;
        I = new Profile(0.0);
    } else {
        // 3. smooth further
        if(n2 > 40) {
            sgSmooth(I2, n2, n2 / 20, n2 / 20, 0, 1);
        }

        // 4. create current profile
        I = new Profile(z2, I2, n2);

        /**5. Normalize profile to match bunch charge as a constant
         * However, only normalize for sufficient beam energy
         */
        thisBunch = this;

        double Qcalc = 0.0;
        double z = zMin;
        dz = (zMax - zMin) / 99;
        for(int i = 1; i < 100; i++) {
            Qcalc += I->get(z, itype_lin);
            z += dz;
        }
        Qcalc *= (dz / (beta * Physics::c));
        I->scale((Q_m > 0.0 ? Q_m : Physics::q_e) / Qcalc);
    }

    dStat |= ds_currentCalculated;

    free(z2);
    free(I2);
    free(b);
    free(z1);
    free(I1);
}

void EnvelopeBunch::cSpaceCharge() {
    Inform msg("cSpaceCharge");

    if(numSlices_m < 2) {
        msg << "EnvelopeBunch::cSpaceCharge called with insufficient slices (" << numMySlices_m << ")" << endl;
        return;
    }

    for(int i = 0; i < numMySlices_m; ++i) {
        Esct[i] = 0.0;
        G[i]    = 0.0;
    }

    if((Q_m <= 0.0) || (I->max() <= 0.0)) {
        msg << "EnvelopeBunch::cSpaceCharge going out" << endl;
        return;
    }

    double Imax = I->max();

    double *xi = vector(numSlices_m);
    int nV = 0;
    double sm = 0.0;
    double A0 = 0.0;

    for(int i = 0; i < numSlices_m; i++) {
        xi[i] = 0.0;
    }
    for(int i = 0; i < numMySlices_m; i++) {
        //if(s[i]->p[SLI_z] > zCat) {
        //nV++;
        A0 = 4.0 * s[i]->p[SLI_x] * s[i]->p[SLI_y];
        sm += A0;
        xi[i+mySliceStartOffset_m] = A0 * (1.0 - s[i]->p[SLI_beta] * s[i]->p[SLI_beta]); // g2
        //}
    }

    //int nVTot = nV;
    //reduce(nVTot, nVTot, OpAddAssign());
    //if(nVTot < 2) {
    //free(xi);
    //return;
    //}
    nV = numMySlices_m;
    int nVTot = numSlices_m;

    reduce(&xi[0], &xi[0] + numSlices_m, &xi[0], OpAddAssign());
    reduce(sm, sm, OpAddAssign());
    A0 = sm / nVTot;

    double dzMin = 5.0 * Physics::c * Q_m / (Imax * numSlices_m);
    for(int localIdx = 0; localIdx < numMySlices_m; localIdx++) {

        double z0 = z_m[localIdx + mySliceStartOffset_m];
        sm = 0.0;

        for(int j = 0; j < numSlices_m; j++) {
            double zj = z_m[j];
            double dz = fabs(zj - z0);
            if((dz > dzMin) && (zj > zCat)) {
                double Aj = xi[j] / (dz * dz);
                double v  = 1.0 - (1.0 / sqrt(1.0 + Aj));
                if(zj > z0) {
                    sm -= v;
                } else {
                    sm += v;
                }
            }
        }

        // longitudinal effect
        double bt = s[localIdx]->p[SLI_beta];
        Esct[localIdx] = (bt < BETA_MIN1 ? 0.0 : bt < BETA_MIN2 ? pow((bt - BETA_MIN1) / (BETA_MIN2 - BETA_MIN1), 2) : 1.0);
        Esct[localIdx] *= Q_m * sm / (Physics::two_pi * Physics::epsilon_0 * A0 * (nVTot - 1));
        G[localIdx] = I->get(z0, itype_lin) / Imax;

        // tweak to compensate for non-relativity
        if(bt < BETA_MIN2) {
            if(s[localIdx]->p[SLI_beta] < BETA_MIN1)
                G[localIdx] = 0.0;
            else
                G[localIdx] *= pow((bt - BETA_MIN1) / (BETA_MIN2 - BETA_MIN1), 2);
        }
    }

    free(xi);
    return;
}

double EnvelopeBunch::moveZ0(double zC) {
    zCat = zC;
    double dz = zC - zHead();
    if(dz > 0.0) {
        for(int i = 0; i < numMySlices_m; i++) {
            s[i]->p[SLI_z] += dz;
        }
        backup(); // save the new state
        *gmsg << "EnvelopeBunch::moveZ0(): bunch moved with " << dz << " m to " << zCat << " m" << endl;
    }

    return dz;
}

void EnvelopeBunch::backup() {
    for(int i = 0; i < numMySlices_m; i++)
        s[i]->backup();
}

double EnvelopeBunch::tReset(double dt) {
    double new_dt = dt;

    if(dt == 0.0) {
        new_dt = t;
        *gmsg << "EnvelopeBunch time reset at z = " << zAvg() << " m with: " << t << " s, new offset: " << t + t_offset << " s";
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
void EnvelopeBunch::derivs(double tc, double Y[], double dYdt[]) {
    double g2    = 1.0 / (1.0 - Y[SLI_beta] * Y[SLI_beta]);
    double g     = sqrt(g2);
    double ecbmg = -1.0 * Physics::q_e * Physics::c * Y[SLI_beta] / (g * Physics::EMASS);
    double c2    = Physics::c * Physics::c;

    /// minimum spot-size due to emittance
    double enxc2 = pow((emtnx0 + emtbx0) * Physics::c / (Y[SLI_beta] * g), 2);
    double enyc2 = pow((emtny0 + emtby0) * Physics::c / (Y[SLI_beta] * g), 2);

    /// transverse space charge
    /// somewhat strange: I expected: c*c*I/(2*Ia) (R. Bakker)
    ///  kpc  = 0.375*c*c*I->max()/Ia;
    double kpc  = 0.5 * Physics::c * Physics::c * I->max() / Physics::Ia;

    // dYdt[SLI_z] = Y[SLI_beta]*c*sqrt(1.0 - (pow(Y[SLI_px0],2) + pow(Y[SLI_py0],2))/pow(c*Y[SLI_beta],2));
    // dYdt[SLI_z] = Y[SLI_beta]*c*cos(Y[SLI_px0]/Y[SLI_beta]/c)*cos(Y[SLI_py0]/Y[SLI_beta]/c);
    dYdt[SLI_z]    = Y[SLI_beta] * Physics::c * cos(sqrt(Y[SLI_px0] * Y[SLI_px0] + Y[SLI_py0] * Y[SLI_py0]) / Y[SLI_beta] / Physics::c);

    //update beta
    //FIXME: -Esl(2)!!!! or -Physics::e0mc?
    dYdt[SLI_beta] = Physics::e0mc * (-Esl(2) + Esct[cS] + Ezw[cS]) / pow(g, 3);

    /// beta * gamma^2 * dbeta/dt
    double bg2dbdt = Y[SLI_beta] * g2 * dYdt[SLI_beta];

    if(solver & sv_radial) {
        dYdt[SLI_x]  = Y[SLI_px];
        dYdt[SLI_px] = (kpc * G[cS] / (Y[SLI_x] * Y[SLI_beta] * pow(g, 3))) +
                       (enxc2 / pow(Y[SLI_x], 3)) - (KRsl(0) * Y[SLI_x]) - (bg2dbdt * Y[SLI_px]);
        dYdt[SLI_y]  = Y[SLI_py];
        dYdt[SLI_py] = (kpc * G[cS] / (Y[SLI_y] * Y[SLI_beta] * pow(g, 3))) +
                       (enyc2 / pow(Y[SLI_y], 3)) - (KRsl(1) * Y[SLI_y]) - (bg2dbdt * Y[SLI_py]);
    } else {
        dYdt[SLI_x]  = Y[SLI_px];
        dYdt[SLI_px] = 0.0;
        dYdt[SLI_y]  = Y[SLI_py];
        dYdt[SLI_py] = 0.0;
    }

    if(solver & sv_offaxis) {
        dYdt[SLI_x0]  = Y[SLI_px0];
        dYdt[SLI_px0] = -KTsl(0) - (bg2dbdt * Y[SLI_px0]) + Physics::e0m * (g * Exw[cS]);
        dYdt[SLI_y0]  = Y[SLI_py0];
        dYdt[SLI_py0] = -KTsl(1) - (bg2dbdt * Y[SLI_py0]) + Physics::e0m * (g * Eyw[cS]);
    } else {
        dYdt[SLI_x0]  = Y[SLI_px0];
        dYdt[SLI_px0] = 0.0;
        dYdt[SLI_y0]  = Y[SLI_py0];
        dYdt[SLI_py0] = 0.0;
    }
}

void EnvelopeBunch::setEx(double emi) {
    emtnx0 = emi;
}

void EnvelopeBunch::setEy(double emi) {
    emtny0 = emi;
}

void EnvelopeBunch::setCharge(double _Q) {
    sign = _Q < 0.0 ? -1 : 1;
    Q_m = abs(_Q);
}

double EnvelopeBunch::AvBField() {
    double bf = 0.0;
    for(int slice = 0; slice < numMySlices_m; slice++) {
        for(int i = 0; i < 3; i++) {
            bf += BF[slice](i);
        }
    }

    reduce(bf, bf, OpAddAssign());
    return bf / numSlices_m;
}

double EnvelopeBunch::AvEField() {
    double ef = 0.0;
    for(int slice = 0; slice < numMySlices_m; slice++) {
        for(int i = 0; i < 3; i++) {
            ef += EF[slice](i);
        }
    }

    reduce(ef, ef, OpAddAssign());
    return ef / numSlices_m;
}

double EnvelopeBunch::Eavg() {
    int nValid = 0;
    double sum = 0.0;
    for(int i = 0; i < numMySlices_m; i++) {
        if((s[i]->p[SLI_z] > zCat) && s[i]->valid) {
            sum += s[i]->gamma();
            nValid++;
        }
    }
    reduce(sum, sum, OpAddAssign());
    reduce(nValid, nValid, OpAddAssign());
    sum /= nValid;
    return (nValid > 0 ? ((Physics::EMASS * Physics::c * Physics::c / Physics::q_e) * (sum - 1.0)) : 0.0);
}

double EnvelopeBunch::get_sPos() {
    // find reference position = centroid?
    double refpos = 0.0;
    size_t count = 0;
    for(int i = 0; i < numMySlices_m; i++) {
        if(s[i]->p[SLI_z] > 0.0) {
            refpos += s[i]->p[SLI_z];
            count++;
        }
    }

    //if(mySliceStartOffset_m <= mid && mySliceEndOffset_m >= mid)
    //refpos = s[mid-mySliceStartOffset_m]->p[SLI_z];

    reduce(refpos, refpos, OpAddAssign());
    reduce(count, count, OpAddAssign());
    return refpos / count;
}

double EnvelopeBunch::zAvg() {
    int nV = 0;
    double sum = 0.0;
    for(int i = 0; i < numMySlices_m; i++) {
        if(s[i]->valid) {
            sum += s[i]->p[SLI_z];
            nV++;
        }
    }

    reduce(nV, nV, OpAddAssign());
    if(nV < 1)
        throw OpalException("EnvelopeBunch", "EnvelopeBunch::zAvg() no valid slices left");

    reduce(sum, sum, OpAddAssign());
    return (sum / nV);
}

double EnvelopeBunch::zTail() {
    double min;

    int i = 0;
    while((i < numMySlices_m) && (!s[i]->valid))
        i++;
    if(i == numMySlices_m)
        throw OpalException("EnvelopeBunch", "EnvelopeBunch::zTail() no valid slices left");
    else
        min = s[i]->p[SLI_z];

    for(i = i + 1; i < numMySlices_m; i++)
        if((s[i]->p[SLI_z] < min) && (s[i]->valid))
            min = s[i]->p[SLI_z];

    reduce(min, min, OpMinAssign());
    return min;
}

double EnvelopeBunch::zHead() {
    double max;

    int i = 0;
    while((i < numMySlices_m) && (s[i]->valid == 0))
        i++;
    if(i == numMySlices_m)
        throw OpalException("EnvelopeBunch", "EnvelopeBunch::zHead() no valid slices left");
    else
        max = s[i]->p[SLI_z];

    for(i = 1; i < numMySlices_m; i++)
        if(s[i]->p[SLI_z] > max) max = s[i]->p[SLI_z];

    reduce(max, max, OpMaxAssign());
    return max;
}

double EnvelopeBunch::getGamma(int i) {
    assert(i < numMySlices_m);
    return s[i]->gamma();
}

double EnvelopeBunch::getBeta(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_beta];
}

void EnvelopeBunch::setZ(int i, double coo) {
    assert(i < numMySlices_m);
    s[i]->p[SLI_z] = coo;
}

double EnvelopeBunch::getZ(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_z];
}

double EnvelopeBunch::getX(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_x];
}

double EnvelopeBunch::getY(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_y];
}

double EnvelopeBunch::getX0(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_x0];
}

double EnvelopeBunch::getY0(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_y0];
}

double EnvelopeBunch::getPx(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_px];
}

double EnvelopeBunch::getPy(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_py];
}

double EnvelopeBunch::getPz(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_beta] * Physics::m_e * s[i]->gamma();
}

double EnvelopeBunch::getPx0(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_px0];
}

double EnvelopeBunch::getPy0(int i) {
    assert(i < numMySlices_m);
    return s[i]->p[SLI_py0];
}

void EnvelopeBunch::updateFields() {
    Inform msg("updateFields()");
    IpplTimings::startTimer(selfFieldTimer_m);

    // Calculate the current profile
    if(Q_m > 0.0) {

        synchronizeSlices();

        calcI();

        // the following assumes space-charges do not change significantly over nSc steps
        cSpaceCharge();

        // accordingly also the wake-fiels are considered static over a single step
        //cWake();

        if(numSlices_m > 40) {
            // smooth Esct to get rid of numerical noise
            //double Esc = 0;
            //sgSmooth(Esct,n,n/20,n/20,0,4);
        }
    } else {
        I = new Profile(0.0);
    }

    IpplTimings::stopTimer(selfFieldTimer_m);
}

void EnvelopeBunch::run(double tStep, double _zCat) {
    Inform msg("tStep");
    static int msgParsed = 0;

    int nok, nbad, i;
    // default accuracy of integration
    double eps = 1.0e-4;
    // integration step-size (tStep by default)
    double dt;
    EnvelopeSlice *sp;

    thisBunch = this;
    zCat     = _zCat;

    int i0 = 0;
    int di = 1;
    int iStart = -1;

    // backup last stage before new execution
    backup();

    dt = tStep;

#ifndef BETEMISSIONMODEL
    if(lastEmittedSlice_m < numSlices_m) {

        dt = tEmission_m;
        size_t nextSlice = numSlices_m - 1 - lastEmittedSlice_m;

        if(mySliceStartOffset_m <= nextSlice && nextSlice <= mySliceEndOffset_m) {
            // emit next slice
            iStart = nextSlice - mySliceStartOffset_m;
            s[iStart]->p[SLI_z] = zCat;
        }

        if(nextSlice > mySliceEndOffset_m)
            iStart = numMySlices_m;

        lastEmittedSlice_m++;
    }
#else
    if(s[0]->p[SLI_z] < zCat) {
        i = numSlices_m - 1;
        while((i > 0) && (s[i]->p[SLI_z] >= zCat)) --i;

        if(i >= 0) {
            iStart = i;
            dt = (zCat - s[i]->p[SLI_z]) / (s[i]->p[SLI_beta] * Physics::c);
            if(dt > tStep) {
                dt = tStep;
                s[iStart]->p[SLI_z] += (Physics::c * dt * s[i]->p[SLI_beta]);
            } else {
                s[iStart]->p[SLI_z]  = zCat;
            }
            for(i = 0; i < iStart; i++) {
                s[i]->p[SLI_z] += (Physics::c * dt * s[i]->p[SLI_beta]);
            }
        }

        lastEmittedSlice_m++;
    }
#endif

    for(i = 0; i < numMySlices_m; i++) {
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
                    rk4(&(sp->p[0]), SLNPAR, t, dt, Gderivs);
                    ode_result = 0;
                } else {
                    ode_result = odeint(&(sp->p[0]), SLNPAR, t, t + dt, epsLocal, 0.1 * dt, 0.0, &nok, &nbad, Gderivs);
                }

                if(ode_result != 0) {
                    // restore the backup
                    sp->restore();
                    epsLocal *= 10.0;
                }
            }

            if(ode_result == 1) {
                // use fixed step integration if dynamic fails
                rk4(&(sp->p[0]), SLNPAR, t, dt, Gderivs);

                if(msgParsed < 2) {
                    msg << "EnvelopeBunch::run() Switched to fixed step RK rountine for solving of DE at slice " << i << endl;
                    msgParsed = 2;
                }
            } else if((epsLocal != eps) && (msgParsed == 0)) {
                msg << "EnvelopeBunch::run() integration accuracy relaxed to " << epsLocal << " for slice " << i << "ONLY FIRST OCCURENCE MARKED!" << endl;
                msgParsed = 1;
            }
        }

        if(s[i]->check()) {
            msg << "Slice " << i << " no longer valid at z = " <<  s[i]->p_old[SLI_z] << " beta = " << s[i]->p_old[SLI_beta] << endl;
            msg << "Slice " << i << " no longer valid at z = " <<  s[i]->p[SLI_z] << " beta = " << s[i]->p[SLI_beta] << endl;
        }
    }
    // mark that slices might not be synchronized (and space charge accordingly)
    dStat &= (!(ds_slicesSynchronized | ds_spaceCharge));

    /// mark calling of this function + update vars
    t += dt;

    /// subtract average orbit for when tracking along the s-axis
    if(solver & sv_s_path) {
        int nV = 0;
        double ga = 0.0, x0a = 0.0, px0a = 0.0, y0a = 0.0, py0a = 0.0;
        double beta, fi_x, fi_y;

        //FIXME: BET calculates only 80 %, OPAL doesn't ?

        for(i = 0; i < numMySlices_m; i++) {
            sp  = s[i];
            if((sp->p[SLI_z] >= zCat) && sp->valid) {
                ++nV;
                ga  += sp->gamma();
                x0a += sp->p[SLI_x0];
                y0a += sp->p[SLI_y0];
                px0a += sp->p[SLI_px0];
                py0a += sp->p[SLI_py0];
            }
        }

        int nVTot = nV;
        reduce(nVTot, nVTot, OpAddAssign());
        if(nVTot > 0) {
            if(nV > 0) {
                reduce(ga, ga, OpAddAssign());
                reduce(x0a, x0a, OpAddAssign());
                reduce(px0a, px0a, OpAddAssign());
                reduce(y0a, y0a, OpAddAssign());
                reduce(py0a, py0a, OpAddAssign());
            }
            ga  = ga / nVTot;
            x0a = x0a / nVTot;
            px0a = px0a / nVTot;
            y0a = y0a / nVTot;
            py0a = py0a / nVTot;
        } else {
            msg << "EnvelopeBunch::run() No valid slices to subtract average" << endl;
        }
        beta = sqrt(1.0 - (1.0 / pow(ga, 2)));
        fi_x = px0a / Physics::c / beta;
        fi_y = py0a / Physics::c / beta;

        dx0 -= x0a;
        dy0 -= y0a;
        dfi_x -= fi_x;
        dfi_y -= fi_y;
        for(i = 0; i < numMySlices_m; i++) {
            sp = s[i];

            sp->p[SLI_x0] -= x0a;
            sp->p[SLI_y0] -= y0a;
            sp->p[SLI_px0] -= px0a;
            sp->p[SLI_py0] -= py0a;
            sp->p[SLI_z] += (sp->p[SLI_x0] * sin(fi_x) + sp->p[SLI_y0] * sin(fi_y));
        }
    }
}

void EnvelopeBunch::initialize(int sli, double Q, double energy, double width, double te, double frac, double current, double center, double bX, double bY, double mX, double mY, double Bz0) {
    createSlices(sli);

    // Set charge Q and energy E0
    this->setCharge(Q);
    this->setEnergy(energy);

    //FIXME: how do I have to set this and width??
    //and why do I have to set zPos of slices??
#ifndef BETEMISSIONMODEL
    center = -1 * te / 2.0;
#else
    //center = -1*te/2.0;
    //te = -te;
    te = 0.0011;
    frac = 0.0;
    center = -0.00055;
#endif
    tEmission_m = te / numSlices_m;
    lastEmittedSlice_m = 0;

    // Shape of the bunch should be rectangular
    int i2 = 0;
    this->setLShape(i2 ? bsGauss : bsRect, center, te, frac);
    this->setTShape(mX, mY, bX, bY, Bz0);

    *gmsg << "tBin = " << tEmission_m << " tEmission = " << te << endl;
    *gmsg << "Q = " << Q << " NSlices =" << sli << " Thermal Energy = " << energy << endl << endl;

    // Set solver method
    // 12: on-axis, radial, default track all)
    // for other modes consult system.C
    this->setSolver(12);
}

// compare BET and OPAL time
// set time of PartBunch to time of EnvelopeBunch
void EnvelopeBunch::actT() {
    setT(this->getT());
}

void EnvelopeBunch::get_bounds(Vector_t &min, Vector_t &max) {
    min[0] = 0.0;
    min[1] = 0.0;
    min[2] = this->zTail();
    max[0] = 0.0;
    max[1] = 0.0;
    max[2] = this->zHead();
}

Inform &EnvelopeBunch::slprint(Inform &os) {
    if(this->getTotalNum() != 0) {  // to suppress Nan's
        os << "* ************** S L B U N C H ***************************************************** " << endl;
        os << "* NSlices= " << this->getTotalNum() << " Qtot= " << Q_m << endl; //" [nC]  Qi= " << abs(qi_m) << " [C]" << endl;
        os << "* Emean= " << get_meanEnergy() * 1e-6 << " [MeV]" << endl;
        os << "* dT= " << this->getdT() << " [s]" << endl;
        os << "* spos= " << this->zAvg() << " [m]" << endl;

        //os << "* rmax= " << rmax_m << " [m]" << endl;
        //os << "* rmin= " << rmin_m << " [m]" << endl;
        //os << "* rms beam size= " << rrms_m << " [m]" << endl;
        //os << "* rms momenta= " << prms_m << " [beta gamma]" << endl;
        //os << "* mean beam size= " << rmean_m << " [m]" << endl;
        //os << "* mean momenta= " << pmean_m << " [beta gamma]" << endl;
        //os << "* rms emmitance= " << eps_m << " (not normalized)" << endl;
        //os << "* rms correlation= " << rprms_m << endl;
        //os << "* tEmission= " << getTEmission() << " [s]" << endl;
        os << "* ********************************************************************************** " << endl;
    }
    return os;
}

