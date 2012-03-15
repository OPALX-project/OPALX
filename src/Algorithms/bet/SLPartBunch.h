#ifdef HAVE_ENVELOPE_SOLVER
#ifndef OPAL_SLPartBunch_HH
#define OPAL_SLPartBunch_HH

// ------------------------------------------------------------------------
// $Rfile: SLPartBunch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class SLPartBunch
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: Andreas Adelmann  and Co. $
//
// ------------------------------------------------------------------------

#include "Algorithms/Particle.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"
#include "Algorithms/bet/EnvelopeSlice.h"
#include "Algorithms/bet/EnvelopeBunch.h"

#include <iosfwd>
#include <vector>
#include "Ippl.h"


class SLPartBunch: public PartBunch {

    const PartData *reference;

public:
    /// Default constructor.
    //  Construct empty bunch.
    SLPartBunch(const PartData *ref);

    /// Conversion.
    SLPartBunch(const std::vector<Particle> &,const PartData *ref);

    SLPartBunch(const SLPartBunch &);
    ~SLPartBunch();

    // needed by OPAL
    long LastSection[1000];  // last em-field section

    // Comes the slice slice
    double                dt[1000];  // eigen-time of the slice
    Vector_t               Z[1000];  // centroid in z of the slice

    Inform &slprint(Inform &os);

    void calcBeamParameters() 
    { 
        //FIXME ???!?
        //INFOMSG("HELLO SL calcBeamParameters" << endl);
    }

    int getLocalNumSlices() { return 1;}
    int getTotalNumSlices() { return getLocalNumSlices();}

    double calcGamma(double b);

    void defineBunch ();

    void iniBetBunch(int sli, double charge, double energy, double width, double frac, double current, double center, double bX, double bY, double mX, double mY, double Bz);

    double getGamma(int i);

    double getBeta(int i);

    void setZ(int i, double zcoo);

    double getZ(int i); // get Z coordinate of ith slice

    double getX(int i); // get X coordinate of ith slice

    double getY(int i); // get Y coordinate of ith slice

    double getX0(int i); // get X0 coordinate of ith slice

    double getY0(int i); // get Y0 coordinate of ith slice

    double getPx(int i);  // get momentum

    double getPy(int i);

    double getPz(int i);

    double getPx0(int i);

    double getPy0(int i);

    void plotR();

    void setKR(Vector_t value, int i);

    void setKT(Vector_t value, int i);

    void setEF(Vector_t value, int i);

    void setBF(Vector_t value, int i);

    int getN();

    void actT();

    void writeBetBunch(const char *fname);

    void writeSlicedBetBunch(const char *fname);

    void writeBetStat(const char *fname);

    void tstep(double cat, double step, size_t nstep);

    void BetOut(FILE* dat, FILE* sli);

public:

    EnvelopeBunch* BetBunch_m; // contains an element bunch from the 'bet' programme

    double wfraction; // fraction of data, that is written in .dat file

};

inline Inform &operator<<(Inform &os, SLPartBunch &p)
{
    return p.slprint(os);
}

#endif // OPAL_SLPartBunch_HH
#endif
