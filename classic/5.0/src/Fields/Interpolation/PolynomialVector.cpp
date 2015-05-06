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



#include <iomanip>
#include <sstream>
#include <math.h>

#include "gsl/gsl_sf_gamma.h"

#include "Utilities/GeneralClassicException.h"

#include "Fields/Interpolation/MMatrix.h"
#include "Fields/Interpolation/MVector.h"

#include "Fields/Interpolation/PolynomialVector.h"


namespace interpolation {

////////// POLYNOMIALVECTOR START ///////////

bool PolynomialVector::_printHeaders=true;

PolynomialVector::PolynomialVector() 
    : _pointDim(0), _polyKeyByPower(),
      _polyKeyByVector(), _polyCoeffs(),
      _polyVector() {
}

PolynomialVector::PolynomialVector (const PolynomialVector& pv) 
    : _pointDim(pv._pointDim), _polyKeyByPower(pv._polyKeyByPower),
      _polyKeyByVector(pv._polyKeyByVector),
      _polyCoeffs(pv._polyCoeffs.num_row(), pv._polyCoeffs.num_col(), 0.),
      _polyVector(pv._polyVector) {
    SetCoefficients(pv._polyCoeffs);
}


PolynomialVector::PolynomialVector(int numberOfInputVariables, MMatrix<double> polynomialCoefficients)
        :  _pointDim(numberOfInputVariables), _polyKeyByPower(), _polyKeyByVector(), _polyCoeffs(polynomialCoefficients)
{
  SetCoefficients(numberOfInputVariables, polynomialCoefficients);
}

PolynomialVector::PolynomialVector(std::vector<PolynomialCoefficient> coefficients)
        :  _pointDim(0), _polyKeyByPower(), _polyKeyByVector(), _polyCoeffs()
{
  SetCoefficients(coefficients);
}

void PolynomialVector::SetCoefficients(int pointDim, MMatrix<double> coeff)
{
  int nPCoeff      = coeff.num_col();
  _pointDim        = pointDim;
  _polyKeyByPower  = std::vector< std::vector<int> >();
  _polyKeyByVector = std::vector< std::vector<int> >();
  _polyVector      = std::vector<PolynomialCoefficient>();

  for(int i=0; i<nPCoeff; i++)
    _polyKeyByPower.push_back(IndexByPower (i, pointDim));
  for(int i=0; i<nPCoeff; i++) 
    _polyKeyByVector.push_back(IndexByVector(i, pointDim));

  for(size_t i=0; i<_polyCoeffs.num_row(); i++) 
    for(int j=0; j<nPCoeff; j++)
      _polyVector.push_back( PolynomialVector::PolynomialCoefficient(_polyKeyByVector[j], i, _polyCoeffs(i+1, j+1) ) );
}

void PolynomialVector::SetCoefficients(std::vector<PolynomialCoefficient> coeff)
{
  _pointDim        = 0;
  _polyKeyByPower  = std::vector< std::vector<int> >();
  _polyKeyByVector = std::vector< std::vector<int> >();
  _polyVector      = coeff;

  int maxPolyOrder = 0;
  int valueDim     = 0;
  for(unsigned int i=0; i<coeff.size(); i++) {
    int polyOrder = coeff[i].InVariables().size();
    for(unsigned int j=0; j<coeff[i].InVariables().size(); j++)
      if(coeff[i].InVariables()[j] > _pointDim) _pointDim = coeff[i].InVariables()[j];
    if(coeff[i].OutVariable()    > valueDim)  valueDim  = coeff[i].OutVariable();
    if(polyOrder > maxPolyOrder) maxPolyOrder = polyOrder;
  }
  _pointDim++; //PointDim indexes from 0
  _polyCoeffs = MMatrix<double>(valueDim+1, NumberOfPolynomialCoefficients(_pointDim, maxPolyOrder+1), 0.);

  for(size_t i=0; i<_polyCoeffs.num_col(); i++)
    _polyKeyByPower.push_back(IndexByPower (i, _pointDim));
  for(size_t i=0; i<_polyCoeffs.num_col(); i++) 
  {
    _polyKeyByVector.push_back(IndexByVector(i, _pointDim));
    for(unsigned int j=0; j<coeff.size(); j++) 
      if(coeff[j].InVariables() == _polyKeyByVector.back()) 
      {
        int dim = coeff[j].OutVariable();
        _polyCoeffs(dim+1,i+1) = coeff[j].Coefficient();
      }
  }
}

void PolynomialVector::SetCoefficients(MMatrix<double> coeff)
{
  for(size_t i=0; i<coeff.num_row() && i<_polyCoeffs.num_row(); i++)
    for(size_t j=0; j<coeff.num_col() && j<_polyCoeffs.num_col(); j++)
      _polyCoeffs(i+1,j+1) = coeff(i+1,j+1);
}

PolynomialVector* PolynomialVector::Recentred(double * point) const
{
    throw(GeneralClassicException("PolynomialVector::Recentred", "Recentred not implemented"));
}

void  PolynomialVector::F(const double*   point,    double* value)          const
{
    MVector<double> pointV(PointDimension(), 1), valueV(ValueDimension(), 1);
    for(unsigned int i=0; i<PointDimension(); i++) pointV(i+1) = point[i];
    F(pointV, valueV);
    for(unsigned int i=0; i<ValueDimension(); i++) value[i]  = valueV(i+1);
}

void  PolynomialVector::F(const MVector<double>& point,   MVector<double>& value) const
{
    MVector<double> polyVector(_polyKeyByVector.size(), 1); 
    MakePolyVector(point, polyVector);
    MVector<double> answer = _polyCoeffs*polyVector;
    for(unsigned int i=0; i<ValueDimension(); i++) value(i+1) = answer(i+1);
}


MVector<double>&  PolynomialVector::MakePolyVector(const MVector<double>& point, MVector<double>& polyVector) const
{
    for(unsigned int i=0; i<_polyKeyByVector.size(); i++)
        for(unsigned int j=0; j<_polyKeyByVector[i].size(); j++)
            polyVector(i+1) *= point( _polyKeyByVector[i][j]+1 );
    return polyVector;
}

double*  PolynomialVector::MakePolyVector(const double* point, double* polyVector) const
{
    for(unsigned int i=0; i<_polyKeyByVector.size(); i++)
    {
        polyVector[i] = 1.;
        for(unsigned int j=0; j<_polyKeyByVector[i].size(); j++)
            polyVector[i] *= point[ _polyKeyByVector[i][j] ];
    }
    return polyVector;
}

double* PolynomialVector::MakeDerivVector(const double* positions, const int* deriv_indices, double* deriv_vec) const {
    size_t pkey_size = _polyKeyByPower.size();
    for (size_t i = 0; i < pkey_size; ++i) {
        deriv_vec[i] = 1.;
        for (int j = 0; j < _pointDim; ++j) {
            int power = _polyKeyByPower[i][j] - deriv_indices[j];
            if (power < 0) {
                deriv_vec[i] = 0.;
            } else {
                // x^(n-m)
                for (int k = 0; k < power; ++k)
                    deriv_vec[i] *= positions[j];
            }
            // n*(n-1)*(n-2)*...*(n-m)
            for (int k = _polyKeyByPower[i][j]; k > 0 && k > power; --k) {
                deriv_vec[i] *= k;
            }
        }
    }
    return deriv_vec;
}

//Turn int index into a std::vector<int> 'a' of length 'd' with values x_1^a_1 x_2^a_2 ... x_d^a_d
// e.g. x_1^4 x_2^3 = {4,3,0,0}
std::vector<int> PolynomialVector::IndexByPower(int index, int nInputVariables)
{
    std::vector<int> powerIndex(nInputVariables, 0);
    std::vector<int> vectorIndex = IndexByVector(index, nInputVariables);
    for(unsigned int i=0; i<vectorIndex.size(); i++)
        powerIndex[vectorIndex[i]]++;
    return powerIndex;
}

//Turn int index into an index for element of polynomial
// e.g. x_1^4 x_2^3 = {1,1,1,1,2,2,2}
std::vector<int> PolynomialVector::IndexByVector(int index, int nInputVariables)
{
    if(index == 0) return std::vector<int>();
    std::vector<int> indices(1,-1);
    nInputVariables--;
    for(int i=0; i<index; i++)
    {
        if     (indices.front() == nInputVariables) { indices = std::vector<int>(indices.size()+1, 0); indices.back()--;}
        if     (indices.back()  == nInputVariables)
        {
            int j=indices.size()-1;
            while(indices[j] == nInputVariables) j--;
            for(int k=indices.size()-1; k>=j; k--) indices[k] = indices[j]+1;
        }
        else    indices.back()++;
    }
    return indices;
}

unsigned int PolynomialVector::NumberOfPolynomialCoefficients(int pointDimension, int order)
{
    int n=0;
    if(order<=0) return 0;
    for(int i=1; i<order; i++)
        n += gsl_sf_choose(pointDimension+i-1, i);
    return n+1;
}

std::ostream& operator<<(std::ostream& out, const PolynomialVector& pv)
{
    if(PolynomialVector::_printHeaders) pv.PrintHeader(out, '.', ' ', 14, true);
    out << '\n' << pv.GetCoefficientsAsMatrix();
    return out;
}

std::ostream& operator << (std::ostream& out, const PolynomialVector::PolynomialCoefficient& coeff) {
    std::vector<int> in_var = coeff.InVariables();
    for (size_t i = 0; i < coeff.InVariables().size(); ++i)
        out << in_var[i] << " ";
    out << "| " << coeff.OutVariable();
    out << " | " << coeff.Coefficient();
    return out;
}

template <class Container >
void PolynomialVector::PrintContainer(std::ostream& out, const Container& container, char T_separator, char str_separator, int length, bool pad_at_start)
{
//  class Container::iterator it;
  std::stringstream strstr1("");
  std::stringstream strstr2("");
  class Container::const_iterator it1 = container.begin();
  class Container::const_iterator it2 = it1; 
  while(it1!=container.end())
  {
    it2++;
    if(it1 != container.end() && it2 != container.end())
         strstr1 << (*it1) << T_separator;
    else strstr1 << (*it1);
    it1++;
  }

  if(int(strstr1.str().size()) > length && length > 0) strstr2 << strstr1.str().substr(1, length);
  else                                                 strstr2 << strstr1.str();

  if(!pad_at_start) out << strstr2.str();
  for(int i=0; i<length - int(strstr2.str().size()); i++) out << str_separator;
  if(pad_at_start) out << strstr2.str();
}

void PolynomialVector::PrintHeader(std::ostream& out, char int_separator, char str_separator, int length, bool pad_at_start) const
{
  if(_polyKeyByPower.size() > 0) PrintContainer< std::vector<int> >(out, _polyKeyByPower[0], int_separator, str_separator, length-1, pad_at_start);
  for(unsigned int i=1; i<_polyKeyByPower.size(); i++)
    PrintContainer<std::vector<int> >(out, _polyKeyByPower[i], int_separator, str_separator, length, pad_at_start);
}

MVector<double> PolynomialVector::Means(std::vector<std::vector<double> > values, std::vector<double> weights)
{
  if(values.size() < 1)    throw(GeneralClassicException("PolynomialVector::Means", "No input values"));
  if(values[0].size() < 1) throw(GeneralClassicException("PolynomialVector::Means", "Dimension < 1"));
  if(weights.size() != values.size()) weights = std::vector<double>(values.size(), 1.);
  int                npoints     = values.size();
  int                dim         = values[0].size();
  MVector<double>    means(dim,0);
  double             totalWeight = 0;
  for(int x=0; x<npoints; x++) totalWeight += weights[x];
  for(int x=0; x<npoints; x++)
    for(int i=0; i<dim; i++)
    {
      means(i+1) += values[x][i]*weights[x]/totalWeight;
    }

  return means;
}

MMatrix<double> PolynomialVector::Covariances(std::vector<std::vector<double> > values, std::vector<double> weights, MVector<double> means)
{
  size_t                 npoints = values.size();
  if(npoints < 1) throw(GeneralClassicException("PolynomialVector::Covariances()", "No input values"));
  size_t                 dim     = values[0].size();
  if(dim < 1)     throw(GeneralClassicException("PolynomialVector::Covariances()", "Dimension < 1"));    
  if(means.num_row() != dim)           means   = Means(values, weights);
  if(weights.size()  != values.size()) weights = std::vector<double>(values.size(), 1.);

  MMatrix<double>    cov(dim,dim,0.);

  double             totalWeight = 0;
  for(size_t x=0; x<npoints; x++)
    totalWeight += weights[x];

  for(size_t x=0; x<npoints; x++)
    for(size_t i=0; i<dim; i++)
      for(size_t j=0; j<dim; j++) //make a sym matrix for speed
        cov(i+1,j+1) += values[x][i]*values[x][j]*weights[x]/totalWeight;

  for(size_t i=0; i<dim; i++)
    for(size_t j=0; j<dim; j++)
      cov(i+1,j+1) -= means(i+1)*means(j+1);

  return cov;
}

void  PolynomialVector::FDeriv(const double*   point, const int* derivative_by_power, double* value) const {
    MVector<double> polyVector(_polyKeyByVector.size(), 1);
    MakeDerivVector(point, derivative_by_power, &polyVector(1));
    MVector<double> answer = _polyCoeffs*polyVector;
    for(unsigned int i=0; i<ValueDimension(); i++)
        value[i] = answer(i+1);
}



PolynomialVector* PolynomialVector::PolynomialLeastSquaresFit(const std::vector< std::vector<double> >& points, const std::vector< std::vector<double> >& values, 
                                                              int polynomialOrder, const std::vector<double>& weights)
{
  //Algorithm: We have F2 = sum_i ( f_k f_l) where f are polynomial terms; FY = sum_i (f_)

  int pointDim = points[0].size();
  int valueDim = values[0].size();
  int nPoints  = points.size();
  int nCoeffs  = NumberOfPolynomialCoefficients(pointDim, polynomialOrder);

  MMatrix<double> Fy(nCoeffs, valueDim, 0.);
  MMatrix<double> F2(nCoeffs, nCoeffs,  0.);

  MMatrix<double> dummy(valueDim, nCoeffs, 0.);
  for(int i=0; i<valueDim; i++) for(int j=0; j<nCoeffs; j++) dummy(i+1, j+1) = 1;
  PolynomialVector* temp = new PolynomialVector(pointDim, dummy);

  std::vector<double>    tempFx(nCoeffs,0); //tempfx = x_i^j
  std::vector<double>    wt    (nPoints,1);
  if(weights.size() > 0) wt = weights; 

  //sum over points, values
  for(int i=0; i<nPoints; i++)
  {
    temp->MakePolyVector(&points[i][0], &tempFx[0]);
    for(int k=0; k<nCoeffs;  k++)
      for(int j=0; j<nCoeffs;  j++)
        F2(j+1, k+1) += tempFx[k]*tempFx[j]*wt[i];
    for(int j=0; j<nCoeffs; j++)
      for(int k=0; k<valueDim; k++)
        Fy(j+1,k+1) += values[i][k]*tempFx[j]*wt[i];
  }

  try { 
    F2.invert();
  } catch(GeneralClassicException exc) {
      delete temp;
      throw(GeneralClassicException(
          "PolynomialVector::PolynomialLeastSquaresFit",
          "Could not find least squares fit for data"));
  }
  MMatrix<double> A = F2 * Fy;
  delete temp;
  temp = new PolynomialVector(pointDim, A.T());
  return temp;
}

PolynomialVector* PolynomialVector::ConstrainedPolynomialLeastSquaresFit(
                  const std::vector< std::vector<double> >&  points, 
                  const std::vector< std::vector<double> >& values, 
                  int polynomialOrder,
                  std::vector< PolynomialVector::PolynomialCoefficient > coeffs,
                  const std::vector<double>& weights) {
  //Algorithm: we want g(x) = old_f(x) + new_f(x),
  //where old_f has polynomial terms in poly, new_f has all the rest
  //Use value - f(point) as input and do LLS forcing old_f(x) polynomial terms to 0
  //Nb: in this constrained case we have to go through the output vector and treat each like a 1D vector
  int pointDim   = points[0].size();
  int valueDim   = values[0].size();
  int nPoints    = points.size();
  if(nPoints<1)                   throw(GeneralClassicException("PolynomialVector::ConstrainedPolynomialLeastSquaresFit(...)", "No data for LLS Fit"));
  if(int(values.size())!=nPoints) throw(GeneralClassicException("PolynomialVector::ConstrainedPolynomialLeastSquaresFit(...)", "Misaligned array in LLS Fit")); 
  int nCoeffsNew = NumberOfPolynomialCoefficients(pointDim, polynomialOrder);
  for (size_t i = 0; i < coeffs.size(); ++i) {
      std::vector<int> in_var = coeffs[i].InVariables();
      if (coeffs[i].OutVariable() >= valueDim)
          throw(GeneralClassicException(
                "PolynomialVector::ConstrainedPolynomialLeastSquaresFit(...)",
                "Coefficient output variable out of range")); 
      if ( int(in_var.size()) >= polynomialOrder) 
          throw(GeneralClassicException(
                "PolynomialVector::ConstrainedPolynomialLeastSquaresFit(...)",
                "Coefficient input variables of too high order")); 
      for (size_t j = 0; j < coeffs[i].InVariables().size(); ++j) {
          if (in_var[j] >= pointDim)
              throw(GeneralClassicException(
                "PolynomialVector::ConstrainedPolynomialLeastSquaresFit(...)",
                "Coefficient input variable of too high dimension")); 
      }
  }

  MMatrix<double> A(valueDim, nCoeffsNew, 0.);

  std::vector<double>    wt    (nPoints,   1);
  if(weights.size() > 0) wt = weights; 

  PolynomialVector  newPolyVector1D(pointDim, MMatrix<double>(1, nCoeffsNew)); //guaranteed to be of right order
  std::vector< std::vector<double> > tempFx(nPoints, std::vector<double>(nCoeffsNew) );
  for(int i=0; i<nPoints; i++)
    newPolyVector1D.MakePolyVector(&points[i][0], &tempFx[i][0]);

  for(int dim=0; dim<valueDim ; dim++)
  {
    std::vector<PolynomialVector::PolynomialCoefficient> oldCoeffs1D;
    for(unsigned int i=0; i<coeffs.size(); i++)
      if(coeffs[i].OutVariable() == dim && int(coeffs[i].InVariables().size()) < polynomialOrder) 
        oldCoeffs1D.push_back(coeffs[i]);
    PolynomialVector oldPolyVector1D(oldCoeffs1D);
    int              nCoeffsOld = oldCoeffs1D.size();
    int              deltaCoeff = nCoeffsNew - nCoeffsOld;
    if(deltaCoeff <= 0) break;


    std::vector<int> needsCalculation; //index of variables that need calculation
    for(int i=0; i<nCoeffsNew; i++) 
    {
      bool exists = true;
      for(unsigned int j=0; j<oldCoeffs1D.size(); j++)
        if(IndexByVector(i, pointDim) == oldCoeffs1D[j].InVariables()) exists = false;
      if(exists) needsCalculation.push_back(i);
    }

    MVector<double> Fy(deltaCoeff, 0);
    MMatrix<double> F2(deltaCoeff, deltaCoeff,  0.);

    for(int i=0; i<nPoints && needsCalculation.size()>0; i++) //optimisation note - this algorithm spends all its time in this loop
    {
      std::vector<double> oldValue(valueDim);
      oldPolyVector1D.F(&points[i][0], &oldValue[0]);
      for(int k=0; k<deltaCoeff;  k++)
      {
        Fy(k+1) += (values[i][dim] - oldValue[dim])*tempFx[i][needsCalculation[k]]*wt[i]; //ynew - yold
        for(int j=0; j<deltaCoeff;  j++) 
          F2(j+1,k+1) += tempFx[i][needsCalculation[k]]*tempFx[i][needsCalculation[j]]*wt[i];
      }
    }
    try{ F2.invert(); }
    catch(GeneralClassicException exc) { 
      throw(GeneralClassicException("PolynomialVector::ConstrainedPolynomialLeastSquaresFit", "Could not find least squares fit for data")); }
    MVector<double> AVec = F2 * Fy;
    for(int i=0; i<deltaCoeff; i++) A(dim+1, needsCalculation[i]+1) = AVec(i+1);
  }

  PolynomialVector tempPoly(coeffs);
  for(unsigned int i=0; i<coeffs.size(); i++)
    for(unsigned int j=0; j<tempPoly._polyKeyByVector.size(); j++)
      if(tempPoly._polyKeyByVector[j] == coeffs[i].InVariables())
      {
        A(coeffs[i].OutVariable()+1,j+1) = coeffs[i].Coefficient();
      }

  return new PolynomialVector(pointDim, A);
}


PolynomialVector* PolynomialVector::Chi2ConstrainedLeastSquaresFit
                                    (const std::vector< std::vector<double> >&  xin, const std::vector< std::vector<double> >& xout, 
                                      int polynomialOrder, std::vector< PolynomialVector::PolynomialCoefficient > coeffs, 
                                      double chi2Start, double discardStep, double* chi2End, double chi2Limit, 
                                      std::vector<double> weights, bool firstIsMean)
{
  int    nParticles   = xin.size();
  int    dimensionOut = xout[0].size();
  if(int(weights.size()) != nParticles) weights = std::vector<double>(xin.size(), 1.);
  std::vector<double> weights_in = std::vector<double>(weights);
  std::vector<double> amps(nParticles);
  MVector<double> means(dimensionOut,0);
  if(!firstIsMean)
     means = Means      (xout, weights);
  else for(int i=0; i<dimensionOut; i++) means(i+1) = xout[0][i];

  MMatrix<double> covInv = Covariances(xout, weights, means);
  try{ covInv.invert(); }
  catch(GeneralClassicException exc) {throw(GeneralClassicException("PolynomialVector::Chi2ConstrainedLeastSquaresFit", "Failed to find least squares fit for data"));}

  double totalWeight = 0.;
  double discard       = chi2Start;
  double chi2          = chi2Limit*2.;
  int    nCoefficients = PolynomialVector::NumberOfPolynomialCoefficients(xin[0].size(), polynomialOrder)*xout[0].size();
  int    nGood         = nParticles;

  if(discard <= 0.) for(int i=0; i<nParticles; i++) //if chi2Start ill-defined, just use maximum chi2 in sample
  {
    MVector<double> xoutVec    (dimensionOut, 0);
    for(int j=0; j<dimensionOut; j++)
      xoutVec(j+1) = xout[i][j] - means(j+1);
    amps[i]      = (xoutVec.T() * covInv * xoutVec)(1);
    if(amps[i] > discard) discard = amps[i];
    totalWeight += weights[i];
  }
  discard /= discardStep; //Set discard to the largest amplitude in the data

  PolynomialVector* pvec = new PolynomialVector(std::vector<PolynomialVector::PolynomialCoefficient>() );
  while( nGood >= nCoefficients && chi2 > chi2Limit)
  {
    chi2               = 0.;
    nGood              = nParticles;
    if(pvec != NULL) delete pvec;

    //optimisation note - algorithm spends all its time doing the CPLS Fit
    try                 { pvec = PolynomialVector::ConstrainedPolynomialLeastSquaresFit(xin, xout, polynomialOrder, coeffs, weights); }
    catch(GeneralClassicException exc) { pvec = NULL; chi2 = chi2Limit*2.;}
    for(int i=0; i<nParticles && pvec!=NULL; i++)
    {
      if(fabs(weights[i]) > 1.e-9) 
      {
        std::vector<double> pout(dimensionOut, 0.);
        pvec->F(&xin[i][0], &pout[0]);
        MVector<double> residualVec(dimensionOut,0);
        for(int j=0; j<dimensionOut; j++) residualVec(j+1) = pout[j] - xout[i][j];
        chi2             += weights[i]*weights[i]/totalWeight/totalWeight * (residualVec.T() * covInv * residualVec)(1); //watch weights handling here
      }
    }
    for(int i=0; i<nParticles; i++)
      if(amps[i] > discard) {totalWeight -= weights[i]; weights[i] = 0.; nGood--;}
    discard /= discardStep;
    std::cout << "ConstrainedPolynomialLeastSquaresFit - chi2: " << chi2 << " cut: " << discard << " cut_step: " << discardStep << " weight: " << totalWeight
              << " pvec: " << pvec << std::endl;
  }
  if(chi2 > chi2Limit || pvec == NULL)
  {
    std::stringstream err;
    err << "PolynomialVector::Chi2ConstrainedLeastSquaresFit failed to converge. Best fit had <chi2>="
                                << chi2 << " compared to limit " << chi2Limit << std::endl;
    throw(GeneralClassicException("PolynomialVector::Chi2ConstrainedLeastSquaresFit(...)", err.str()));
  }
  if(chi2End != NULL) *chi2End = discard; //save discard at end for future use
  return pvec;
}

void PolynomialVector::PolynomialCoefficient::SpaceTransform(std::vector<int> space_in, std::vector<int> space_out)
{
  std::map<int, int> mapping; //probably optimise this
  for(unsigned int i=0; i<space_out.size(); i++)
    for(unsigned int j=0; j<space_in.size(); j++)
      if(space_out[i] == space_in[j]) mapping[j] = i; //mapping[space_in_index] returns space_out_index

  std::vector<int> in_variables(_inVarByVec.size());
  for(unsigned int con=0; con<in_variables.size(); con++)
  {
    if(mapping.find(in_variables[con]) != mapping.end()) in_variables[con] = mapping[in_variables[con]];
    else                                                throw(GeneralClassicException("PolynomialVector::PolynomialCoefficient::SpaceTransform", "Input variable not found in space transform"));
  }

  if(mapping.find(_outVar) != mapping.end()) _outVar = mapping[_outVar];
  else                                      throw(GeneralClassicException(
                                                                  "PolynomialVector::PolynomialCoefficient::SpaceTransform",
                                                                  "Output variable not found in space transform"
));

  _inVarByVec = in_variables;
}
////////// POLYNOMIALVECTOR END ////////////

}


