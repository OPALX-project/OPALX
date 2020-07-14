//
// Class: RealDiracMatrix
//   Real Dirac matrix class.
//   They're ordered after the paper of Dr. C. Baumgarten: "Use of real Dirac matrices in two-dimensional coupled linear optics".
//   The diagonalizing method is based on the paper "Geometrical method of decoupling" (2012) of Dr. C. Baumgarten.
//
// Copyright (c) 2014, 2020 Christian Baumgarten, Paul Scherrer Institut, Villigen PSI, Switzerland
//                          Matthias Frey, ETH Zürich
// All rights reserved
//
// Implemented as part of the Semester thesis by Matthias Frey
// "Matched Distributions in Cyclotrons"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef RDM_H
#define RDM_H

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

class RealDiracMatrix
{
public:
    /// Sparse matrix type definition
    typedef boost::numeric::ublas::compressed_matrix<
        double,
        boost::numeric::ublas::row_major
    > sparse_matrix_t;

    /// Dense matrix type definition
    typedef boost::numeric::ublas::matrix<double> matrix_t;
    /// Dense vector type definition
    typedef boost::numeric::ublas::vector<double> vector_t;

    /// Default constructor (sets only NumOfRDMs and DimOfRDMs)
    RealDiracMatrix();

    /*!
        * @param i specifying the matrix (has to be in the range from 0 to 15)
        * @returns the i-th Real Dirac matrix
        */
    sparse_matrix_t getRDM(short);

    /*!
     * Decomposes a real-valued 4x4 matrix into a linear combination
     * and returns a vector containing the coefficients
     * @param M an arbitrary real-valued 4x4 matrix
     */
    vector_t decompose(const matrix_t&);

    /*!
     * Takes a vector of coefficients, evaluates the linear combination of
     * RDMs with these coefficients and returns a 4x4 matrix
     *
     * @param coeffs is a vector of coefficients (at most length NumOfRDMs)
     */
    matrix_t combine(const vector_t&);

    /*!
     * Brings a 4x4 symplex matrix into Hamilton form and
     * computes the transformation matrix and its inverse
     *
     * @param Ms is a 4x4 symplex matrix
     * @param R is the 4x4 transformation matrix (gets computed)
     * @param invR is the 4x4 inverse transformation matrix (gets computed)
     */
    void diagonalize(matrix_t&, sparse_matrix_t&, sparse_matrix_t&);

    /*!
     * Diagonalizes (in-place) the 4x4 sigma matrix. This algorithm
     * is Implemented according to https://arxiv.org/abs/1205.3601
     *
     * @param sigma is the 4x4 sigma matrix
     * @returns the inverse of the total transformation
     */
    sparse_matrix_t diagonalize(matrix_t&);

    /*!
     * @param M 4x4 real-valued matrix
     * @returns the symplex part of a 4x4 real-valued matrix
     */
    matrix_t symplex(const matrix_t&);

    /*!
     * @param M 4x4 real-valued matrix
     * @returns the cosymplex part of a 4x4 real-valued matrix
     */
    matrix_t cosymplex(const matrix_t&);

    /// The number of real Dirac matrices
    short NumOfRDMs;
    /// The matrix dimension (4x4)
    short DimOfRDMs;


private:
    ///
    /*!
     * Applies a rotation to the matrix M by a given angle
     * @param M is the matrix to be transformed
     * @param i is the i-th RDM used for transformation
     * @param phi is the angle of rotation
     * @param Rtot is a reference to the current transformation matrix
     * @param invRtot is a reference to the inverse of the current transformation matrix
     */
    void transform(matrix_t&, short, double, sparse_matrix_t&, sparse_matrix_t&);

    void transform(short, double, sparse_matrix_t&, sparse_matrix_t&);

    void update(matrix_t& sigma, short i, sparse_matrix_t& R, sparse_matrix_t& invR);
};

#endif