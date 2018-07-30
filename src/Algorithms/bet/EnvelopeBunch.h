#ifndef _ENVELOPE_BUNCH_H
#define _ENVELOPE_BUNCH_H

#include "Algorithms/bet/profile.h"
#include "Algorithms/bet/EnvelopeSlice.h"
#include "Utilities/OpalException.h"
#include "Algorithms/PartBunch.h"
#include <Physics/Physics.h>

#include <assert.h>
#include <vector>

#include <memory>

enum EnvelopeBunchParameter {
    sp_beta,      /// normalized velocity (total) [-]
    sp_gamma,     /// Lorenz factor
    sp_z,         /// slice position [m]
    sp_I,         /// current [A]
    sp_Rx,        /// beam size x [m]
    sp_Ry,        /// beam size y [m]
    sp_Px,        /// beam divergence x
    sp_Py,        /// beam divergence y
    sp_Pz,        /// beam divergence z
    sp_x0,        /// position centroid x [m]
    sp_y0,        /// position centroid y [m]
    sp_px0,       /// angular deflection centriod x
    sp_py0        /// angular deflection centroid y
};

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

enum EnvelopeBunchShape {
    bsRect,
    bsGauss
};


/**
 * @brief core of the envelope tracker based on Rene Bakkers BET
 * implementation
 */
class EnvelopeBunch : public PartBunch {

public:
    /// Default constructor
    EnvelopeBunch(const PartData *ref);

    /// Conversion
    EnvelopeBunch(const std::vector<OpalParticle> &,
                  const PartData *ref);

    /// Copy constructor (implement if needed)
    EnvelopeBunch(const EnvelopeBunch &) = delete;

    virtual ~EnvelopeBunch();

    /// create and initialize local num slices
    void createBunch();

    /// distributes nSlice amongst processors and initializes slices
    void distributeSlices(int nSlice = 101);

    bool isValid_m;

    /// current profile of bunch (fit)
    std::unique_ptr<Profile> currentProfile_m;

    IpplTimings::TimerRef calcITimer_m;
    IpplTimings::TimerRef spaceChargeTimer_m;

    /// initialize an envelope bunch (interface for Distribution class)
    using PartBunch::initialize;
    void initialize(int sli, double charge, double energy, double width, double te, double frac,
                    double current, double center, double bX, double bY, double mX, double mY, double Bz, int nbin);

    /// check if solver includes radial
    bool isRadial() { return solver & sv_radial; }

    /// check if solver includes off-axis tracking
    bool isOffaxis() { return solver & sv_offaxis; }

    /// calculates envelope statistics
    void calcBeamParameters();

    void computeSpaceCharge();

    /**
     * @brief performs a time-step for all active slices (also handles
     * emission)
     *
     * @param tStep dt for timestep
     * @param zCat position of cathode (default: 0.0)
     */
    void timeStep(double tStep, double zCat = 0.0);

    /**
     * @brief helper function to calculate derivatives need in RH equation
     *
     * @param tc time
     * @param Y[] in/out array of slice parameters
     * @param dYdt[] derivation wrt to time of slice parameters
     */
    void derivs(double tc, double Y[], double dYdt[]);

    /// calculate the average energy of the bunch
    double Eavg();
    /// calculate <z> [m]
    double zAvg();
    /// calculate tail of bunch [m]
    double zTail();
    /// calculate the head of the bunch [m]
    double zHead();
    /// read time-stamp of bunch
    double time() { return t; }

    /// set emittance X
    void setEx(double emi) { emtnx0 = emi; }
    /// set emittance Y
    void setEy(double emi) { emtny0 = emi; }
    /// set the DE solver flag
    void setSolverParameter(int s) { solver = s; }
    // set particle energy of bunch in [eV] and optional the correlated energy spread [eV/m]
    void setEnergy(double, double = 0.0);

    void setExternalFields(int i, Vector_t EF, Vector_t BF, Vector_t KR, Vector_t KT) {
        this->EF[i] = EF;
        this->BF[i] = BF;
        this->KR[i] = KR;
        this->KT[i] = KT;
    }

    /// return reference position
    double get_sPos();
    /// returns average magnetic field
    double AvBField();
    /// returns average electric field
    double AvEField();
    /// returns the number of local slices
    int getLocalNum() { return numMySlices_m; }
    /// returns the total number of slices
    int getTotalNum() { return numSlices_m; }
    /// returns the current time of the bunch
    double getT() { return t; }
    /// returns the mean energy
    double get_meanKineticEnergy() { return Eavg(); }
    /// returns the energy spread
    double get_dEdt() { return dEdt_m; }
    /// returns vector with rms position
    Vector_t sigmax() { return sigmax_m; }
    /// returns vector with rms momenta
    Vector_t sigmap() { return sigmap_m; }
    /// returns vector with emittance
    Vector_t emtn() { return emtn_m; }
    /// returns vector with normalized emittance
    Vector_t get_norm_emit() { return norm_emtn_m; }
    /// returns vector with the max spatial extends of the bunch
    Vector_t maxX() { return maxX_m; }
    /// returns vector with the min spatial extends of the bunch
    Vector_t minX() { return minX_m; }
    /// returns vector with the max momentum of the bunch
    Vector_t maxP() { return maxP_m; }
    /// returns vector with the min momentum of the bunch
    Vector_t minP() { return minP_m; }
    /// returns charge per slice
    double getChargePerParticle() { return Q_m / numSlices_m; }
    /// returns RMS x,y,z
    Vector_t get_rrms() { return sigmax_m; }
    /// returns RMSP x,y,z
    Vector_t get_prms() { return sigmap_m; }

    size_t mySliceStartOffset() { return mySliceStartOffset_m; }
    size_t mySliceEndOffset() { return mySliceEndOffset_m; }
    size_t numMySlices() { return numMySlices_m; }

    Inform &slprint(Inform &os);

    /// set the charge of the bunch
    void setCharge(double _Q) {
        sign = _Q < 0.0 ? -1 : 1;
        Q_m = std::abs(_Q);
    }

    /// returns gamma of slice i
    double getGamma(int i) {
        assert(i < numMySlices_m);
        return s[i]->computeGamma();
    }

    /// returns beta of slice i
    double getBeta(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_beta];
    }
    void setBeta(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_beta] = val;
    }

    /// set Z coordinate of slice i
    void setZ(int i, double coo) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_z] = coo;
    }

    /// returns Z coordinate of slice i
    double getZ(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_z];
    }

    /// returns X coordinate of slice i
    double getX(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_x];
    }
    void setX(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_x] = val;
    }

    /// returns Y coordinate of slice i
    double getY(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_y];
    }
    void setY(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_y] = val;
    }

    /// returns X coordinate of the centroid of slice i
    double getX0(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_x0];
    }
    void setX0(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_x0] = val;
    }

    /// returns Y coordinate of the centroid of slice i
    double getY0(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_y0];
    }
    void setY0(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_y0] = val;
    }

    /// returns X momenta of slice i
    double getPx(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_px];
    }
    void setPx(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_px] = val;
    }

    /// returns Y momenta of slice i
    double getPy(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_py];
    }
    void setPy(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_py] = val;
    }

    /// returns Z momenta of slice i
    double getPz(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_beta] * Physics::m_e * s[i]->computeGamma();
    }

    /// returns angular deflection centroid in x of slice i
    double getPx0(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_px0];
    }
    void setPx0(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_px0] = val;
    }

    /// returns angular deflection centroid in y of slice i
    double getPy0(int i) {
        assert(i < numMySlices_m);
        return s[i]->p[SLI_py0];
    }
    void setPy0(int i, double val) {
        assert(i < numMySlices_m);
        s[i]->p[SLI_py0] = val;
    }

    /// returns bounds of envelope bunch
    void get_bounds(Vector_t &min, Vector_t &max) {
        min[0] = 0.0;
        min[1] = 0.0;
        min[2] = this->zTail();
        max[0] = 0.0;
        max[1] = 0.0;
        max[2] = this->zHead();
    }


private:
    /// number of total slices in bunch
    int numSlices_m;
    /// number of my slices in bunch
    int numMySlices_m;
    /// first global slice on this processor
    size_t mySliceStartOffset_m;
    /// last global slice on this processor
    size_t mySliceEndOffset_m;
    /// synchronized z positions for parallel tracker
    std::vector<double> z_m;
    /// synchronized betas for parallel tracker
    std::vector<double> b_m;
    /// bins for emission
    std::vector< std::vector<int> > bins_m;
    /// emission bin width
    double hbin_m;
    /// number of bins for emission
    int nebin_m;
    /// first bin on processor containing slices
    int firstBinWithValue_m;
    /// number of active slices
    int activeSlices_m;

    /// see enum SolverParameters
    int solver;
    /// see enum DataStatus
    int dStat;
    /// local time in bunch [s]
    double t;
    /// accumulated time offset by tReset function
    double t_offset;
    /// intrinsic normalized emittance of slice [m rad]
    double emtnx0, emtny0;
    /// intrinsic normalized emittance Bush effect [m rad]
    double emtbx0, emtby0;
    /// offset of the coordinate system when tracking along the s-axis [m]
    double dx0, dy0;
    /// rotation of coordinate system when tracking along the s-axis [rad]
    double dfi_x, dfi_y;
    /// magnetic field on cathode [T]
    double Bz0;
    /// total bunch charge [C]
    double Q_m;
    /// average current on creation of bunch (see setLshape)
    double I0avg;
    /// electric field
    Vector_t Esl;
    /// magnetic field
    Vector_t Bsl;
    /// radial focussing term beam
    Vector_t KRsl;
    /// transverse kick of beam
    Vector_t KTsl;
    /// define value of radial kick for each slice
    std::unique_ptr<Vector_t[]> KR;
    /// define value of transversal kick for each slice
    std::unique_ptr<Vector_t[]> KT;
    /// external E fields
    std::unique_ptr<Vector_t[]> EF;
    /// external B fields
    std::unique_ptr<Vector_t[]> BF;
    /// array of slices
    std::vector< std::shared_ptr<EnvelopeSlice> > s;
    /// gives the sign of charge Q
    int sign;
    /// current Slice set in run() & cSpaceCharge() and used in derivs() & zcsI()
    int cS;
    /// cathode position
    double zCat;
    /// transverse wake field x
    std::vector<double> Exw;
    /// transverse wake field y
    std::vector<double> Eyw;
    /// longitudinal wake field
    std::vector<double> Ezw;
    /// Longitudinal Space-charge field
    std::vector<double> Esct;
    /// Transverse Space-charge term: Eq.(9)
    std::vector<double> G;

    int nValid_m;
    double z0_m;
    double emission_time_step_;
    size_t lastEmittedBin_m;
    double E_m;
    double dEdt_m;
    double Einc_m;
    double tau_m;
    double I_m;
    double Irms_m;
    double Rx_m;
    double Ry_m;
    double RxMax_m;
    double RyMax_m;
    double RxMin_m;
    double RyMin_m;
    double Px_m;
    double Py_m;
    double x0_m;
    double y0_m;
    double x0Max_m;
    double y0Max_m;
    double x0Min_m;
    double y0Min_m;
    double dx0_m;
    double dy0_m;
    double dfi_x_m;
    double dfi_y_m;
    double Ez_m;
    double Bz_m;
    Vector_t maxX_m;
    Vector_t minX_m;
    Vector_t maxP_m;
    Vector_t minP_m;
    Vector_t sigmax_m;
    Vector_t sigmap_m;
    Vector_t emtn_m;
    Vector_t norm_emtn_m;

    // used in derivs() to remove call to zHead/Tail() containing MPI
    // collectives causing problems if slices are unevenly distributed across
    // processors.
    double curZHead_m;
    double curZTail_m;

    /**
     * @brief synchronize z position and betas of all slices (needed in calcI
     * and space charge calculation)
     */
    void synchronizeSlices();

    /**
     * @brief run statistics on slices
     *
     * @param sp parameter to run statistics on
     * @param xAvg average
     * @param xMax max
     * @param xMin min
     * @param rms rms
     * @param nValid number of valid slices
     */
    void runStats(EnvelopeBunchParameter sp, double *xAvg, double *xMax, double *xMin, double *rms, int *nValid);

    /**
     * @brief calculate bunch emittance
     *
     * @param emtnx normalized emittance x
     * @param emtny normalized emittance y
     * @param emtx emittance x
     * @param emty emittance y
     * @param nValid number of valid slices
     */
    void calcEmittance(double *emtnx, double *emtny, double *emtx, double *emty, int *nValid);

    /**
     * @brief calculate the energy chirp and uncorrelated energy spread
     *
     * @param g0 average gamma
     * @param dgdt chirp
     * @param gInc incoherent energy spread
     * @param nValid number of valid slices
     */
    void calcEnergyChirp(double *g0, double *dgdt, double *gInc, int *nValid);

    /**
     * @brief set longitudinal shape of bunch (initial distribution)
     *
     * @param shape of bunch (currently rectangular or Gauss)
     * @param z0 center of the bunch [m]
     * @param w length of the bunch [m]
     * @param frac fraction of Gauss (length) used
     */
    void setBinnedLShape(EnvelopeBunchShape shape, double z0, double w, double frac);

    /**
     * @brief set transverse shape of bunch (initial distribution)
     *
     * @param enx normalized emittance x [m rad]
     * @param eny normalized emittance y [m rad]
     * @param rx radius x [m]
     * @param ry radius y [m]
     * @param b0 Bz0 [T]
     */
    void setTShape(double enx, double eny, double rx, double ry, double b0);

    /**
     * @brief set transverse offset of bunch
     *
     * @param x0 coordinate [m]
     * @param px0 divergence [rad]
     * @param y0 coordinate [m]
     * @param py0 divergence [rad]
     */
    void setTOffset(double x0, double px0, double y0, double py0);

    /// move the complete bunch forward such that the head of the bunch matches the cahtode position
    double moveZ0(double zC);

    /// backup slice values
    void backup() {
        for (auto & slice : s)
            slice->backup();
    }


    /// reset time of bunch (returns the offset applied)
    /// time difference (0.0 - auto-sync)
    double tReset(double dt = 0.0);

    /**
      * @briefcal culates space-charge fields
      * Calculate longitudinal-, and transverse space-charge fields for bunch
      * Output is stored in the global class arrays Esct and G, respectively.
      */
    void cSpaceCharge();

    /// calculates the current current distribution
    void calcI();
};

inline Inform &operator<<(Inform &os, EnvelopeBunch &p) {
    return p.slprint(os);
}
#endif

