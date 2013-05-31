// ------------------------------------------------------------------------
// $RCSfile: StaticFixedPoint.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class StaticFixedPoint:
//   Fixed point for class Vps.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/StaticFixedPoint.h"
#include "Algebra/Matrix.h"
#include "Algebra/LUMatrix.h"
#include "Algebra/TpsMonomial.h"
#include "Algebra/Vector.h"
#include "Algebra/VpsInvMap.h"
#include <cmath>


// Tolerance for fixed point search.
namespace {
    const double tol = 1.0e-10;
}


// Class StaticFixedPoint
// ------------------------------------------------------------------------

StaticFixedPoint::StaticFixedPoint():
    fixedPoint(),
    fixedPointMap()
{}


StaticFixedPoint::StaticFixedPoint(const StaticFixedPoint &fp):
    fixedPoint(fp.fixedPoint),
    fixedPointMap(fp.fixedPointMap)
{}


StaticFixedPoint::StaticFixedPoint(const VpsInvMap<double> &map):
    fixedPoint(VpsInvMap<double>::identity(map.getDimension())),
    fixedPointMap(map.getDimension()) {
    const int nFree = 4;
    const int nDim = map.getDimension();
    LUMatrix<double> lu;
    Vector<double> fixPoint(nDim, 0.0);

    // Find fixed point;
    // tempMap is the map around the fixed point.
    VpsInvMap<double> tempMap(map);

    while(true) {
        // Get system of equations for fixed point.
        Matrix<double> A = tempMap.linearTerms();
        Vector<double> Error = tempMap.constantTerm();
        double error = 0.0;

        for(int i = 0; i < nFree; i++) {
            A(i, i) -= 1.0;
            if(std::abs(Error(i)) > error) error = std::abs(Error(i));
        }

        for(int i = nFree; i < nDim; i++) {
            for(int j = 0; j < nDim; j++) A(i, j) = A(j, i) = 0.0;
            A(i, i) = 1.0;
            Error(i) = 0.0;
        }

        // Correction for fixed point.
        lu = LUMatrix<double>(A);
        lu.backSubstitute(Error);
        fixPoint -= Error;

        // Test for convergence.
        if(error < tol) break;

        // Build map around fixed point found so far.
        tempMap = map.substitute(fixPoint + fixedPoint) - fixPoint;
    }

    // Fixed point map for orders 1 and higher in delta (dispersion).
    for(int limit = 1; limit <= map.getTopOrder(); limit++) {
        // Remove pure delta terms.
        VpsInvMap<double> Q = tempMap.substitute(fixedPoint) -
                              fixedPoint.substitute(fixedPointMap);
        TpsMonomial delta(nDim);
        delta[nFree+1] = limit;
        Vector<double> Error(nDim, 0.0);
        for(int i = 0; i < nFree; i++) {
            fixedPointMap[i] += Q[i].filter(limit, limit);
            Error(i) = Q[i].getCoefficient(delta);
            fixedPointMap[i].setCoefficient(delta, 0.0);
        }

        lu.backSubstitute(Error);

        for(int i = 0; i < nFree; i += 2) {
            // Effect on orbit.
            fixedPoint[i].setCoefficient(delta, - Error(i));
            fixedPoint[i+1].setCoefficient(delta, - Error(i + 1));

            // Orbit terms acting on path length difference.
            TpsMonomial powers(nDim);
            powers[nFree+1] = limit - 1;
            powers[i] = 1;
            fixedPoint[nFree].setCoefficient(powers, double(limit) * Error(i + 1));
            powers[i] = 0;
            powers[i+1] = 1;
            fixedPoint[nFree].setCoefficient(powers, - double(limit) * Error(i));
            powers[i+1] = 0;
        }

        // Path length term in delta^limit.
        fixedPoint[nFree].
        setCoefficient(delta, Q[nFree].getCoefficient(delta));
    }

    // Fill in path length terms and identity for delta(p)/p.
    fixedPointMap[nFree+1] = Tps<double>::makeVariable(nDim, nFree + 1);
    fixedPointMap[nFree] =
        fixedPoint.inverse()[nFree].substitute(tempMap).substitute(fixedPoint);

    // Combine fixed point and dispersion.
    fixedPoint += fixPoint;
}


StaticFixedPoint::~StaticFixedPoint()
{}


const VpsInvMap<double> &StaticFixedPoint::getFixedPoint() const {
    return fixedPoint;
}


const VpsInvMap<double> &StaticFixedPoint::getFixedPointMap() const {
    return fixedPointMap;
}
