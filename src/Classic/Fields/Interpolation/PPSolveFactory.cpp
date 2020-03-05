/* 
 *  Copyright (c) 2015, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer. 
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to 
 *     endorse or promote products derived from this software without specific 
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>

#include "Utility/Inform.h" // ippl

#include "Utilities/GeneralClassicException.h"

#include "Fields/Interpolation/SolveFactory.h"
#include "Fields/Interpolation/PolynomialPatch.h"
#include "Fields/Interpolation/PPSolveFactory.h"

extern Inform* gmsg;

namespace interpolation {

typedef PolynomialCoefficient Coeff;

namespace {
    std::vector<double> getDeltaPos(Mesh* mesh) {
        Mesh::Iterator it = mesh->begin();
        std::vector<double> pos1 = it.getPosition();
        for (size_t i = 0; i < it.getState().size(); ++i)
            it[i]++;
        std::vector<double> pos2 = it.getPosition();
        for (size_t i = 0; i < pos2.size(); ++i)
            pos2[i] -= pos1[i];
        return pos2;
    }
}

PPSolveFactory::PPSolveFactory(Mesh* points,
                               std::vector<std::vector<double> > values,
                               int polyPatchOrder,
                               int smoothingOrder) 
  : polyPatchOrder_m(polyPatchOrder),
    smoothingOrder_m(smoothingOrder),
    polyDim_m(0),
    polyMesh_m(),
    points_m(points),
    values_m(values),
    polynomials_m(),
    thisPoints_m(),
    thisValues_m(),
    derivPoints_m(),
    derivValues_m(),
    derivIterator_m(),
    derivIndices_m() {
    if (points == NULL) {
        throw GeneralClassicException(
                          "PPSolveFactory::Solve",
                          "NULL Mesh input to Solve");
    }
    if (points->end().toInteger() != int(values.size())) {
        std::stringstream ss;
        ss << "Mismatch between Mesh and values in Solve. Mesh size was "
           << points->end().toInteger() << " but field values size was "
           << values.size();
        throw GeneralClassicException(
                          "PPSolveFactory::Solve",
                          ss.str());
    }
    if (polyPatchOrder < 1) {
        throw GeneralClassicException(
                "PPSolveFactory::Solve",
                "Polynomial order must be at least 1.");
    }
    if (smoothingOrder < polyPatchOrder) {
        throw GeneralClassicException(
                "PPSolveFactory::Solve",
                "Polynomial order must be <= smoothing order in Solve");
    }
    polyMesh_m = points->dual();
    if (polyMesh_m == NULL)
        throw GeneralClassicException(
                          "PPSolveFactory::Solve",
                          "Dual of Mesh was NULL");
    polyDim_m = values[0].size();
}

void PPSolveFactory::getPoints() {
    int posDim = points_m->getPositionDimension();
    std::vector< std::vector<int> > dataPoints =
                          getNearbyPointsSquares(posDim, -1, polyPatchOrder_m);
    thisPoints_m = std::vector< std::vector<double> >(dataPoints.size());
    std::vector<double> deltaPos = getDeltaPos(points_m);
    // now iterate over the index_by_power_list elements - offset; each of
    // these makes a point in the polynomial fit
    for (size_t i = 0; i < dataPoints.size(); ++i) {
        thisPoints_m[i] = std::vector<double>(posDim);
        for (int j = 0; j < posDim; ++j)
            thisPoints_m[i][j] = (0.5-dataPoints[i][j])*deltaPos[j];
     }
}

void PPSolveFactory::getValues(Mesh::Iterator it) {
    thisValues_m = std::vector< std::vector<double> >();
    int posDim = it.getState().size();
    std::vector< std::vector<int> > dataPoints =
                          getNearbyPointsSquares(posDim, -1, polyPatchOrder_m);
    size_t dataPointsSize = dataPoints.size();
    Mesh::Iterator end = points_m->end()-1;
    // now iterate over the indexByPowerList elements - offset; each of
    // these makes a point in the polynomial fit
    for (size_t i = 0; i < dataPointsSize; ++i) {
        std::vector<int> itState = it.getState();
        for (int j = 0; j < posDim; ++j)
            itState[j]++;
        Mesh::Iterator itCurrent = Mesh::Iterator(itState, points_m);
        bool outOfBounds = false; // element is off the edge of the mesh
        for (int j = 0; j < posDim; ++j) {
            itCurrent[j] -= dataPoints[i][j];
            outOfBounds = outOfBounds ||
                            itCurrent[j] < 1 ||
                            itCurrent[j] > end[j];
        }
        if (outOfBounds) { // if off the edge, then just constrain to zero
            thisValues_m.push_back(values_m[it.toInteger()]);
        } else { // else fit using values
            thisValues_m.push_back(values_m[itCurrent.toInteger()]);
        }
     }
}

/*
 *  calculate a list of Mesh::iterators, each corresponding to the offset of the
 *  point at which the derivative should be calculated relative to the central
 *  polynomial. Goes hand-in-hand with derivIndices.
 *  for e.g. 2d mesh with polyPatchOrder 2, derivatives are always taken from the
 *  4th row/column:
 *  P o o x
 *  o o o x
 *  o o o x
 *  x x x x
 *  and derivatives used are for smoothOrder 4
 *  P . . - + - +
 *  . . . - + - +
 *  . . . - + - +
 *  | | | \
 *  + + +   +
 *  | | |     \
 *  + + +       +
 *  where each "line" means take a derivative; diagonal line means take a double
 *  derivative (just like in the values, we take all derivatives up to
 *  d_1^n d_2^n for an "nth order" PolynomialSquareVector)
 * 
 *  Note that # marks the polynomial being fitted - derivatives are taken a few
 *  mesh points away from the polynomial (which may not be the ideal). Also note 
 *  that the fit is one-sided; polynomials to the left of # are tasked with
 *  matching points to the right, possibly including #. # knows nothing about
 *  them. We might get a different solution if we reversed the order of fitting,
 *  which is a bit ugly.
 */
void PPSolveFactory::getDerivPoints() {
    int posDim = points_m->getPositionDimension();
    std::vector< std::vector<int> > indices = getNearbyPointsSquares(posDim,
                                                  polyPatchOrder_m,
                                                  smoothingOrder_m);
    std::cerr << "PPSolveFactory::getDerivPoints" << std::endl;
    // derivIterator_m points at the polynomial on which the derivative should
    // be operated (relative offset to centre) - the "x" in comment above
    derivIterator_m = std::vector< Mesh::Iterator >();
    // derivative indexes the order of the derivative - the "+" in the comment
    // above
    derivIndices_m = std::vector< std::vector<int> >();
    derivPoints_m = std::vector< std::vector<double> >();
    for (auto deriv: indices) {
        std::vector<int> derivPointState(deriv);
        for (auto& i: derivPointState) {
            if (i > polyPatchOrder_m) {
                i = polyPatchOrder_m;
            }
        }
        Mesh::Iterator it(derivPointState, polyMesh_m);
        derivPoints_m.push_back(it.getPosition());
        derivIterator_m.push_back(it);
        
        for (auto& i: deriv) {
            if (i > polyPatchOrder_m) {
                i -= polyPatchOrder_m;
            } else {
                i = 0;
            }
        }
        derivIndices_m.push_back(deriv);

    }
    getDerivPolyVec();
    std::cerr << "PPSolveFactory::getDerivPoints 2" << std::endl;
    for (size_t i = 0; i < derivIndices_m.size(); ++i) {
        std::cerr << "    ";
        for (auto n: indices[i]) {
            std::cerr << n << " ";
        }
        std::cerr << " * ";
        for (auto n: derivIndices_m[i]) {
            std::cerr << n << " ";
        }
        std::cerr << " * ";
        for (auto n: derivIterator_m[i].getState()) {
            std::cerr << n << " ";
        }
        std::cerr << " * ";
        for (auto x: derivPoints_m[i]) {
            std::cerr << x << " ";
        }
        std::cerr << std::endl;
    }
    return;
}
    
void PPSolveFactory::getDerivPolyVec() {
    int posDim = points_m->getPositionDimension();
    int nPolyCoeffs = SquarePolynomialVector::NumberOfPolynomialCoefficients
                                                    (posDim, smoothingOrder_m);
    std::vector<double> deltaPos = getDeltaPos(points_m);

    for (size_t i = 0; i < derivPoints_m.size(); ++i) {
        derivPolyVec_m.push_back(MVector<double>(nPolyCoeffs, 1.));
        // Product[x^{p_j-q_j}*p_j!/(p_j-q_j)!]
        for (int j = 0; j < nPolyCoeffs; ++j) {
            std::vector<int> pVec =
                               SquarePolynomialVector::IndexByPower(j, posDim);
            std::vector<int> qVec = derivIndices_m[i];
            for (int l = 0; l < posDim; ++l) {
                int pow = pVec[l]-qVec[l];
                if (pow < 0) {
                    derivPolyVec_m.back()(j+1) = 0.;
                    break;
                }
                derivPolyVec_m.back()(j+1) *= gsl_sf_fact(pVec[l])/gsl_sf_fact(pow)*gsl_sf_pow_int(-deltaPos[l]/2., pow);
            }
        }
    }
}

// Aim is to extract derivatives in the region of "it" and place them into 
// derivValues_m, where the derivValues_m correspond to the derivatives at
// each point in derivPoints_m
void PPSolveFactory::getDerivs(Mesh::Iterator it) {
    // deltaOrder corresponds to the number of derivatives required
    int deltaOrder = smoothingOrder_m - polyPatchOrder_m;
    if (deltaOrder <= 0) {
        return; // no derivatives need to be found
    }

    // derivValues_m is filled with the calc'd derivatives
    derivValues_m = std::vector< std::vector<double> >(derivIterator_m.size());

    std::cerr << "  PPSolveFactory::getDerivs" << std::endl;
    // now get the outer layer of points and put derivatives in for each point
    for (size_t i = 0; i < derivIterator_m.size(); ++i) {
        Mesh::Iterator derivPoint = derivIterator_m[i];
        derivPoint.addState(it);
        if (derivPoint.isOutOfBounds()) {
            // off the edge of the grid - use boundary condition
            std::cerr << " OOB" << std::endl;
            derivValues_m[i] = std::vector<double>(polyDim_m, 0.);
            continue;
        }
        SquarePolynomialVector* poly = polynomials_m[derivPoint.toInteger()];
        if (poly == NULL) {
            // polynomial has not been found yet - this is an error condition
            // we should only use a boundary condition when we are off the mesh;
            // otherwise, we should have calculated a polynomial
            std::cerr << "While working on    " << it << std::endl;
            std::cerr << "  Failed with point " << derivPoint << std::endl;
            std::cerr << "  DerivIterator was " << derivIterator_m[i] << std::endl;
            throw(GeneralClassicException("PPSolveFactory::getDerivs",
                                    "Internal error in polynomial calculation"));
        }
        // polynomial has been found; calculate derivatives
        std::cerr << " Not NULL 0" << std::endl;
        MMatrix<double> coeffs = poly->GetCoefficientsAsMatrix();
        std::cerr << " Not NULL " << derivPolyVec_m.size() << " ** "
                  << coeffs.num_row() << " " << coeffs.num_col() << " ** "
                  << derivPolyVec_m[i].num_row() << std::endl;
        MVector<double> values = coeffs*derivPolyVec_m[i];
        std::cerr << " Not NULL 2" << std::endl;
        derivValues_m[i] = std::vector<double>(polyDim_m);
        std::cerr << " Not NULL 3" << std::endl;
        for(int j = 0; j < polyDim_m; ++j) {
            std::cerr << " Not NULL 4" << std::endl;
            derivValues_m[i][j] = values(j+1);
        }
    }
}

PolynomialPatch* PPSolveFactory::solve() {
    Mesh::Iterator begin = points_m->begin();
    Mesh::Iterator end = points_m->end()-1;
    int meshSize = polyMesh_m->end().toInteger();
    polynomials_m = std::vector<SquarePolynomialVector*>(meshSize, NULL);
    // get the list of points that are needed to make a given poly vector
    getPoints();
    getDerivPoints();
    SolveFactory solver(smoothingOrder_m,
                        polyMesh_m->getPositionDimension(),
                        polyDim_m,
                        thisPoints_m,
                        derivPoints_m,
                        derivIndices_m);
    int total = (polyMesh_m->end()-1).toInteger();
    double oldPercentage = 0.;
    for (Mesh::Iterator it = polyMesh_m->end()-1;
         it >= polyMesh_m->begin();
         --it) {
        double newPercentage = (total-it.toInteger())/double(total)*100.;
        if (newPercentage - oldPercentage > 10.) {
            *gmsg << "    Done " << newPercentage << " %" << endl;
            oldPercentage = newPercentage;
        }
        std::cerr << "PPSolveFactory::solve at " << it << std::endl;
        // find the set of values and derivatives for this point in the mesh
        getValues(it);
        getDerivs(it);
        // The polynomial is found using simultaneous equation solve
        polynomials_m[it.toInteger()] = solver.PolynomialSolve
                                                  (thisValues_m, derivValues_m);
        std::cerr << "  Done " << polynomials_m[it.toInteger()] << " " << it << std::endl;
    }
    return new PolynomialPatch(polyMesh_m, points_m, polynomials_m);
}

void PPSolveFactory::nearbyPointsRecursive(
                               std::vector<int> check,
                               size_t checkIndex,
                               size_t polyPower,
                               std::vector<std::vector<int> >& nearbyPoints) {
    check[checkIndex] = polyPower;
    nearbyPoints.push_back(check);
    if (checkIndex+1 == check.size())
        return;
    for (int i = 1; i < int(polyPower); ++i)
        nearbyPointsRecursive(check, checkIndex+1, i, nearbyPoints);
    for (int i = 0; i <= int(polyPower); ++i) {
        check[checkIndex] = i;
        nearbyPointsRecursive(check, checkIndex+1, polyPower, nearbyPoints);
    }
}

std::vector<std::vector<int> > PPSolveFactory::getNearbyPointsSquares(
                                                      int pointDim,
                                                      int polyOrderLower,
                                                      int polyOrderUpper) {
    if (pointDim < 1)
        throw(GeneralClassicException("PPSolveFactory::getNearbyPointsSquares",
                                      "Point dimension must be > 0"));
    if (polyOrderLower > polyOrderUpper)
        throw(GeneralClassicException("PPSolveFactory::getNearbyPointsSquares",
                   "Polynomial lower bound must be <= polynomial upper bound"));
    // polyOrder -1 (or less) means no terms
    if (polyOrderUpper == polyOrderLower || polyOrderUpper < 0)
        return std::vector<std::vector<int> >();
    // polyOrder 0 means constant term
    std::vector<std::vector<int> > nearbyPoints(1, std::vector<int>(pointDim, 0));
    // polyOrder > 0 means corresponding poly terms (1 is linear, 2 is quadratic...)
    for (size_t thisPolyOrder = 0;
         thisPolyOrder < size_t(polyOrderUpper);
         ++thisPolyOrder) 
        nearbyPointsRecursive(nearbyPoints[0], 0, thisPolyOrder+1, nearbyPoints);
    if (polyOrderLower < 0)
        return nearbyPoints; // no terms in lower bound
    int nPointsLower = 1;
    for (int dim = 0; dim < pointDim; ++dim)
        nPointsLower *= polyOrderLower+1;
    nearbyPoints.erase(nearbyPoints.begin(), nearbyPoints.begin()+nPointsLower);
    return nearbyPoints;
}
}
