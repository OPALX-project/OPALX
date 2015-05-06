// MAUS WARNING: THIS IS LEGACY CODE.
// Copyright 2010-2011 Chris Rogers
// 
// This file is a part of G4MICE
//
// G4MICE is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// G4MICE is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with G4MICE in the doc folder.  If not, see 
// <http://www.gnu.org/licenses/>.

#ifndef PolynomialVector_hh
#define PolynomialVector_hh 1

#include "MMatrix.h"
#include "MVector.h"

#include <map>

/// PolynomialVector, an arbitrary order polynomial vector class.

/// PolynomialVector describes a vector of multivariate polynomials \f$y_i = a_0 + Sum (a_j x^j)\f$
/// i.e. maps a vector \f$\vec{x}\f$ onto a vector \f$\vec{y}\f$ with \f$\vec{y} = a_0 + sum( a_{j_1j_2...j_n} x_1^{j1} x_2^{j2} ... x_n^{jn}  )\f$.
/// Also algorithms to map a single index J into \f$j_1...j_n\f$.\n
/// \n
/// PolynomialVector represents the polynomial coefficients as a matrix of doubles so that\n
/// \f$ \vec{y} = \mathbf{A} \vec{w} \f$\n
/// where \f$\vec{w}\f$ is a vector with elements given by \n
/// \f$ w_i = x_1^{i_1} x_2^{i_2} ... x_n^{i_n} \f$ \n
/// So the index \f$ i \f$ is actually itself a vector. The vectorisation of \f$ i \f$ is handled by IndexByPower and IndexByVector methods.
/// IndexByPower gives an index like:\n
/// \f$ w_i = x_1^{i_1} x_2^{i_2} ... x_n^{i_n} \f$ \n
/// While IndexByVector gives an index like:\n
/// \f$ w_i = x_{i_1}x_{i_2} \ldots x_{i_n} \f$ \n
/// \n
/// Polynomial Vector includes several functions to do least squares fits.
/// Here we find a polynomial of arbitrary dimension and order from a set of points by the method of least squares.
/// Note that the problem must be overspecified, so the number of points must be >= number of polynomial coefficients.
/// A few interesting variations on this theme... have a look!\n
///\n

class PolynomialVector
{
public:

    class PolynomialCoefficient;
    /// Default constructor (leaves everything empty, call SetCoefficients to set up properly)
    PolynomialVector();
    /// Copy constructor
    PolynomialVector (const PolynomialVector& pv); 
    /// Construct polynomial vector passing polynomial coefficients as a matrix.
    PolynomialVector (int pointDimension, MMatrix<double> polynomialCoefficients); 
    /// Construct polynomial vector passing polynomial coefficients as a list of PolynomialCoefficient objects.

    /// Any coefficients missing from the vector are set to 0. Maximum order of the polynomial is given by the 
    /// maximum order of the coefficients in the vector.
    PolynomialVector (std::vector<PolynomialCoefficient> coefficients); 
    /// Destructor - no memory allocated so doesn't do anything
    ~PolynomialVector() {;}

    /// Reinitialise the polynomial vector with new point (x) dimension and coefficients
    void             SetCoefficients(int pointDim, MMatrix<double> coeff);
    /// Reinitialise the polynomial vector with new point (x) dimension and coefficients.

    /// Any coefficients missing from the vector are set to 0. Maximum order of the polynomial is given by the 
    /// maximum order of the coefficients in the vector.
    void             SetCoefficients(std::vector<PolynomialCoefficient> coeff);
    /// Set the coefficients.

    /// This method can't be used to change point dimension or 
    /// value dimension - any range mismatch will be ignored.
    void             SetCoefficients(MMatrix<double> coeff);
    /// Return the coefficients as a matrix of doubles.
    MMatrix<double>                    GetCoefficientsAsMatrix() const {return _polyCoeffs;}
    /// Return the coefficients as a vector of PolynomialCoefficients, including 0 values.
    std::vector<PolynomialCoefficient> GetCoefficientsAsVector() const {return _polyVector;}

    /// Fill value with \f$ y_i \f$ at some set of \f$ x_i \f$ (point).

    /// point should be array of length PointDimension().
    /// value should be array of length ValueDimension().
    void  F    (const double*   point,    double* value)                  const;

    /** Calculate derivative of the field at a given point
     */
    void  FDeriv(const double*   point, const int* derivative_by_power, double* value) const;


    /// Fill value with \f$ y_i \f$ at some set of \f$ x_i \f$ (point).

    /// point should be vector of length PointDimension().
    /// value should be vector of length ValueDimension().
    /// Note that there is no bound checking here.
    void  F    (const MVector<double>& point,   MVector<double>& value) const;
    /// Length of the input point (x) vector.
    unsigned int PointDimension()         const {return _pointDim;}
    /// Length of the output value (y) vector.
    unsigned int ValueDimension()         const {return _polyCoeffs.num_row();}
    /// Index of highest power - 0 is const, 1 is linear, 2 is quadratic etc...
    unsigned int PolynomialOrder()        const {if(_polyKeyByVector.size()>0) return _polyKeyByVector.back().size()+1; else return 0;}
    /// Polymorphic copy constructor.

    /// This is a special copy constructor for inheritance structures
    PolynomialVector* Clone()             const {return new PolynomialVector(*this);}
    /// Return a copy, centred on point.
    PolynomialVector* Recentred(double* point) const;
    /// Return operator $\f R = P(Q) \f$ - NOT IMPLEMENTED
    PolynomialVector  Chain(PolynomialVector Q) const;
    /// Return inverse of the polynomial truncated at order n (in general an infinite series that does not converge, so beware!) - NOT IMPLEMENTED
    PolynomialVector  Inverse(int max_order) const;

    /// Make a vector like \f$(c, x, x^2, x^3...)\f$.
    /// polyVector should be of size NumberOfPolynomialCoefficients().
    /// could be static but faster as member function (use lookup table for _polyKey).
    MVector<double>& MakePolyVector(const MVector<double>& point, MVector<double>& polyVector) const;
    /// Make a vector like \f$(c, x, x^2, x^3...)\f$.
    /// PolyVector should be of size NumberOfPolynomialCoefficients().
    /// Could be static but faster as member function (use lookup table for _polyKey).
    double*          MakePolyVector(const double* point, double* polyVector) const;

    /** Make a vector like \f$d^\vec{p}(c, x, x^2, x^3...)/d\vec{x}^\vec{p}\f$.
     *  - positions: array of size PointDimension() containing the position at
     *    which the vector should be calculated, \f$\vec{x}\f$
     *  - deriv_indices: array in the "index by power" format defining the
     *    derivative index \vec{p}.
     *  - deriv_vec: array that will hold the return value; should be
     *    initialised to size at least NumberOfPolynomialCoefficients()
     */
    double* MakeDerivVector(const double* positions, const int* deriv_indices, double* deriv_vec) const;

    /// Transforms from a 1d index of polynomial coefficients to an nd index.
    /// This is slow - you should use it to build a lookup table.
    /// For polynomial term \f$x_1^i x_2^j ... x_d^n\f$ index like [i][j] ... [n]
    /// e.g. \f$x_1^4 x_2^3 = x_1^4 x_2^3 x_3^0 x_4^0\f$ = {4,3,0,0}.
    static std::vector<int> IndexByPower (int index, int nInputVariables);
    /// Transforms from a 1d index of polynomial coefficients to an nd index.
    /// This is slow - you should use it to build a lookup table.
    /// For polynomial term \f$x_i x_j...x_n\f$ index like [i][j] ... [n]
    /// e.g. \f$x_1^4 x_2^3\f$ = {1,1,1,1,2,2,2}
    static std::vector<int> IndexByVector(int index, int nInputVariables);
    /// Returns the number of coefficients required for an arbitrary dimension, order polynomial
    /// e.g. \f$a_0 + a_1 x + a_2 y + a_3 z + a_4 x^2 + a_5 xy + a_6 y^2 + a_7 xz + a_8 yz + a_9 z^2\f$
    /// => NumberOfPolynomialCoefficients(3,2) = 9 (indexing starts from 0).
    static unsigned int     NumberOfPolynomialCoefficients(int pointDimension, int order);
    unsigned int            NumberOfPolynomialCoefficients() {return NumberOfPolynomialCoefficients(_pointDim, PolynomialOrder());}

    /// Write out the PolynomialVector (effectively just polyCoeffs).
    friend std::ostream& operator<<(std::ostream&,  const PolynomialVector&);

    /// Control the formatting of the polynomial vector. If PrintHeaders is true, then every time I write
    /// a PolynomialVector it will write the header also (default).
    static void PrintHeaders(bool willPrintHeaders) {_printHeaders = willPrintHeaders;}
    /// Write out the header (PolynomialByPower index for each element).
    void PrintHeader(std::ostream& out, char int_separator='.', char str_separator=' ', int length=14, bool pad_at_start=true) const;

    /// Return a PolynomialVector found from a set of points by method of least squares.

    /// Finds the polynomial pol that minimises sum_i ( w_i (\vec{y_i} - Pol(\vec{x_i}))^2
    /// y_i = values[i], x_i = points[i], w_i = weightFunction(x_i) OR weights[i], polynomialOrder is the order of pol
    /// Nb if w_i not defined, I assume that w_i = 1 i.e. unweighted
    /// points and values should be the same length; points[i] should all be the same length; values[i] should all be the same length
    /// weight function should be a vector map that has PointDimension of points and returns a weight of dimension 1
    static PolynomialVector* PolynomialLeastSquaresFit(const std::vector< std::vector<double> >&  points, const std::vector< std::vector<double> >& values, 
                                                       int polynomialOrder, const std::vector<double>& weights=std::vector<double>());
    /// Find a polynomial least squares fit given that I know certain coefficients already, (stored in coeffs)
    /// (e.g. I know the polynomial to 2nd order, I want to extend to 3rd order using a least squares)
    static PolynomialVector* ConstrainedPolynomialLeastSquaresFit(const std::vector< std::vector<double> >&  points, const std::vector< std::vector<double> >& values, 
                                                       int polynomialOrder, std::vector< PolynomialCoefficient > coeffs,
                                                       const std::vector<double>& weights=std::vector<double>());
    /// Find a polynomial least squares fit given that I know certain coefficients already
    /// After finding the fit, calculate "chi2" for each particle. Compare x_out with polynomial fit of x_in, x_pol, 
    /// using delta_chi2 = \vec(x_out - x_pol).T() V^{-1} \vec(x_out - x_pol). If sum(delta_chi2)/totalWeight > chi2Limit,
    /// calculate chi2_out = \vec(x_out).T() V^{-1} \vec(x_out) and discard particles with chi2_out > chi2Start*(discardStep^n)
    /// repeat until sum(delta_chi2)/totalWeight <= chi2Limit, each time incrementing n by 1
    /// if chi2Start < 0, start with maximum chi2; if chi2End != NULL, put the final value in chi2End
    static PolynomialVector* Chi2ConstrainedLeastSquaresFit
                                                      (const std::vector< std::vector<double> >&  xin, const std::vector< std::vector<double> >& xout, 
                                                       int polynomialOrder, std::vector< PolynomialVector::PolynomialCoefficient > coeffs, 
                                                       double chi2Start, double discardStep, double* chi2End, double chi2Limit, 
                                                       std::vector<double> weights, bool firstIsMean=false);
/*
    static PolynomialVector* PolynomialSolve(
             int polynomialOrder,
             const std::vector< std::vector<double> >& positions,
             const std::vector< std::vector<double> >& values,
             const std::vector< std::vector<double> > &deriv_positions,
             const std::vector< std::vector<double> >& deriv_values,
             const std::vector< std::vector<int> >& deriv_indices);
*/
    /// Now some utility functions
    /// Should probably go in a "utility function" namespace, that does not currently exist.

    /// Return the mean of some set of values, using some set of weights
    static MVector<double> Means      (std::vector<std::vector<double> > values, std::vector<double> weights);
    /// Return the covariance matrix of some set of values, using some set of weights and some mean
    /// If length of weights is not equal to length of values use weights = 1.
    /// If length of mean is not equal to length of first value recalculate mean
    static MMatrix<double> Covariances(std::vector<std::vector<double> > values, std::vector<double> weights, MVector<double> mean);
    /// Return chi2 of residuals of some set of output variables compared with the polynomial vector of some set of input variables 
    /// Here Chi2 Avg is defined as Sum( v^T M^{-1} v ) where v is the vector of residuals and M is covariance matrix of v
    double GetAvgChi2OfDifference(std::vector< std::vector<double> > in, std::vector< std::vector<double> > out);
    /// Should probably go in a "utility function" namespace, that does not currently exist.
    /// Print a sequence: with some string that separates elements of index and some character that pads the total vector
    /// length is the total length of the output, set to < 1 to ignore (means no padding)
    /// pad_at_start determines whether to pad at start (i.e. right align) or end (i.e. left align)
    template < class Container >
    static void PrintContainer(std::ostream& out, const Container& container, char T_separator, char str_separator, int length, bool pad_at_start);
    /// Return true if a and b are identical; a must have an iterator object dereferenced by *it
    /// must have increment prefix operator defined and with != operator defined by the object pointed to by it
    template <class TEMP_CLASS, class TEMP_ITER>
    static bool IterableEquality(const TEMP_CLASS& a, const TEMP_CLASS& b, TEMP_ITER a_begin, TEMP_ITER a_end, TEMP_ITER b_begin, TEMP_ITER b_end);

    /// PolynomialCoefficient sub-class - useful, could be extended if needed (e.g. become an iterator etc)
    class PolynomialCoefficient
    {
      public:
        /// inVariablesByVector is x power indexed like e.g. \f$x_1^4 x_2^3 =\f$ {1,1,1,1,2,2,2}
        PolynomialCoefficient(std::vector<int> inVariablesByVector, int outVariable, double coefficient) 
            : _inVarByVec(inVariablesByVector), _outVar(outVariable), _coefficient(coefficient) {}
        PolynomialCoefficient(const PolynomialCoefficient& pc)
            : _inVarByVec(pc._inVarByVec), _outVar(pc._outVar), _coefficient(pc._coefficient) {}
        /// Accessors
        /// Pass a parameter to *set*
        std::vector<int> InVariables() const {return _inVarByVec;}
        int              OutVariable() const {return _outVar;}
        double           Coefficient() const {return _coefficient;}

        std::vector<int> InVariables(std::vector<int> inVar ) {_inVarByVec  = inVar;  return _inVarByVec;}
        int              OutVariable(int              outVar) {_outVar      = outVar; return _outVar;}
        double           Coefficient(double           coeff ) {_coefficient = coeff;  return _coefficient;}

        /// transform coefficient from subspace space_in to subspace space_out, both subspaces of some larger space
        /// if any of coeff variables is not in space_out OR not in space_in, leave this coeff untouched and Squeal
        /// so for coeff({1,2},0,1.1), coeff.space_transform({0,2,3,5}, {4,7,1,2,3,0}) would return coeff({3,4},5,1.1)
        void             SpaceTransform(std::vector<int> space_in, std::vector<int> space_out);
        /// if var is in inVariables return true
        bool             HasInVariable(int var) const {for(unsigned int i=0; i<_inVarByVec.size(); i++) if(_inVarByVec[i] == var) return true; return false;}



      private:
        std::vector<int> _inVarByVec;
        int              _outVar;
        double           _coefficient;
    };


private:
    int                                _pointDim;
    std::vector< std::vector<int> >    _polyKeyByPower; //std::vector<int>[i_1] = j_1
    std::vector< std::vector<int> >    _polyKeyByVector; //std::vector<int>[j_1] = i_1
    MMatrix<double>                    _polyCoeffs;
    std::vector<PolynomialCoefficient> _polyVector;
    static bool                        _printHeaders;
};

std::ostream& operator<<(std::ostream&, const PolynomialVector&);
std::ostream& operator << (std::ostream&, const PolynomialVector::PolynomialCoefficient&);

//template function

template <class TEMP_CLASS, class TEMP_ITER>
bool PolynomialVector::IterableEquality(const TEMP_CLASS& a, const TEMP_CLASS& b, TEMP_ITER a_begin, TEMP_ITER a_end, TEMP_ITER b_begin, TEMP_ITER b_end)
{
  TEMP_ITER a_it=a_begin;
  TEMP_ITER b_it=b_begin;
  while(a_it!=a_end && b_it!=b_end) {
    if( *a_it != *b_it ) return false;
    ++a_it;
    ++b_it;
  }
  if ( a_it != a_end || b_it != b_end ) return false;
  return true;
}



#endif




