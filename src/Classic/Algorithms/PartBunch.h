#ifndef OPAL_PartBunch_HH
#define OPAL_PartBunch_HH

// ------------------------------------------------------------------------
// $RCSfile: PartBunch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class PartBunch
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: Andreas Adelmann  and Co. $
//
// ------------------------------------------------------------------------

#include "Ippl.h"
#include "Algorithms/PBunchDefs.h"
#include "Algorithms/Particle.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartBinsCyc.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Quaternion.h"
#include "Utilities/SwitcherError.h"
#include "Physics/Physics.h"

#include <iosfwd>
#include <vector>


class Distribution;
class LossDataSink;
class FieldSolver;
class ListElem;

template <class T, int, int> class FMatrix;
template <class T, int> class FVector;

// Class PartBunch.
// ------------------------------------------------------------------------
/// Particle Bunch.
//  A representation of a particle bunch as a vector of particles.

// class PartBunch: public std::vector<Particle>, public IpplParticleBase< ParticleSpatialLayout<double, 3> > {
class PartBunch: public IpplParticleBase< ParticleSpatialLayout<double, 3> > {

public:
    /// Default constructor.
    //  Construct empty bunch.
    PartBunch(const PartData *ref);

    /// Conversion.
    PartBunch(const std::vector<Particle> &, const PartData *ref);

    PartBunch(const PartBunch &);
    ~PartBunch();

    void runTests();

    double getRho(int x, int y, int z);

    // MATTHIAS CHECK
    void calcLineDensity(unsigned int nBins, std::vector<double> &lineDensity, std::pair<double, double> &meshInfo);

    /*

      Mesh and Field Layout related functions

    */

    // MATTHIAS CHECK
    const Mesh_t &getMesh() const;

    // MATTHIAS CHECK
    Mesh_t &getMesh();

    // MATTHIAS CHECK
    FieldLayout_t &getFieldLayout();

    void setBCAllPeriodic();
    void setBCAllOpen();

    void setBCForDCBeam();

    /*
      Compatibility function push_back

    */
    
    Inform &print(Inform &os);
    
    std::pair<Vector_t, Vector_t> PartBunch::getEExtrema();
    
    void   computeSelfFields();

    /** /brief used for self fields with binned distribution */
    void   computeSelfFields(int b);

    void computeSelfFields_cycl(double gamma);
    void computeSelfFields_cycl(int b);
    
    // MATTHIAS CHECK
    void resetInterpolationCache(bool clearCache = false);
    
    void swap(unsigned int i, unsigned int j);
    
    
    /// scalar potential
    Field_t rho_m;

    /// vector field on the grid
    VField_t  eg_m;

private:
    
    void updateDomainLength();
    
    void updateFields(const Vector_t& hr, const Vector_t& origin);
    
    /// resize mesh to geometry specified
    void resizeMesh();

    

    bool interpolationCacheSet_m;

    ParticleAttrib<CacheDataCIC<double, 3U> > interpolationCache_m;

    

    PartBunch &operator=(const PartBunch &) = delete;
};



inline
double PartBunch::getRho(int x, int y, int z) {
    return rho_m[x][y][z].get();
}

inline
const Mesh_t &PartBunch::getMesh() const {
    return getLayout().getLayout().getMesh();
}

inline
Mesh_t &PartBunch::getMesh() {
    return getLayout().getLayout().getMesh();
}

inline
FieldLayout_t &PartBunch::getFieldLayout() {
    return dynamic_cast<FieldLayout_t &>(getLayout().getLayout().getFieldLayout());
}

inline
void PartBunch::resetInterpolationCache(bool clearCache) {
    interpolationCacheSet_m = false;
    if(clearCache) {
        interpolationCache_m.destroy(interpolationCache_m.size(), 0, true);
    }
}

#endif // OPAL_PartBunch_HH