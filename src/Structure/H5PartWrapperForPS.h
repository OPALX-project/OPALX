#ifndef OPAL_H5PARTWRAPPERFORPS_H
#define OPAL_H5PARTWRAPPERFORPS_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  H5PartWrapperForPS: a class that manages all calls to H5Part for the ParallelSlice tracker
*/

#include "Structure/H5PartWrapper.h"

#include "Algorithms/PBunchDefs.h"
#include "Utilities/OpalException.h"

#include "H5hut.h"

class EnvelopeBunch;

class H5PartWrapperForPS: public H5PartWrapper {
public:
    H5PartWrapperForPS(const std::string &fileName, h5_int32_t flags = H5_O_WRONLY);
    H5PartWrapperForPS(const std::string &fileName, int restartStep, std::string sourceFile, h5_int32_t flags = H5_O_RDWR);
    virtual ~H5PartWrapperForPS();

    virtual void readHeader();
    virtual void readStep(PartBunch&, h5_ssize_t firstParticle, h5_ssize_t lastParticle);

    virtual void writeHeader();
    virtual void writeStep(PartBunch&, const std::map<std::string, double> &additionalStepAttributes);

    virtual bool predecessorIsSameFlavour() const;

    void stashPhaseSpaceEnvelope(EnvelopeBunch &bunch,
                                 Vector_t FDext[],
                                 double sposHead,
                                 double sposRef,
                                 double sposTail);
    void dumpStashedPhaseSpaceEnvelope();

private:
    void readStepHeader(PartBunch&);
    void readStepData(PartBunch&, h5_ssize_t, h5_ssize_t);

    void writeStepHeader(PartBunch&, const std::map<std::string, double> &);
    void writeStepData(PartBunch&);

    std::vector< Vektor<h5_float64_t, 3> > stash_RefPartR;
    std::vector< Vektor<h5_float64_t, 3> > stash_RefPartP;
    std::vector< Vektor<h5_float64_t, 3> > stash_centroid;
    std::vector< Vektor<h5_float64_t, 3> > stash_geomvareps;
    std::vector< Vektor<h5_float64_t, 3> > stash_xsigma;
    std::vector< Vektor<h5_float64_t, 3> > stash_psigma;
    std::vector< Vektor<h5_float64_t, 3> > stash_vareps;
    std::vector< Vektor<h5_float64_t, 3> > stash_rmin;
    std::vector< Vektor<h5_float64_t, 3> > stash_rmax;
    std::vector< Vektor<h5_float64_t, 3> > stash_maxP;
    std::vector< Vektor<h5_float64_t, 3> > stash_minP;
    std::vector< Vektor<h5_float64_t, 3> > stash_Bhead;
    std::vector< Vektor<h5_float64_t, 3> > stash_Ehead;
    std::vector< Vektor<h5_float64_t, 3> > stash_Bref;
    std::vector< Vektor<h5_float64_t, 3> > stash_Eref;
    std::vector< Vektor<h5_float64_t, 3> > stash_Btail;
    std::vector< Vektor<h5_float64_t, 3> > stash_Etail;
    std::vector< h5_float64_t > stash_actPos;
    std::vector< h5_float64_t > stash_t;
    std::vector< h5_float64_t > stash_meanEnergy;
    std::vector< h5_float64_t > stash_mass;
    std::vector< h5_float64_t > stash_charge;
    std::vector< h5_float64_t > stash_sposHead;
    std::vector< h5_float64_t > stash_sposRef;
    std::vector< h5_float64_t > stash_sposTail;
    std::vector< size_t > stash_nLoc;
    std::vector< size_t > stash_nTot;
};

inline
bool H5PartWrapperForPS::predecessorIsSameFlavour() const {
    return (predecessorOPALFlavour_m == "opal-t");
}

#endif //OPAL_H5PARTWRAPPERFORPS_H