#ifndef MAD_NormalForm_HH
#define MAD_NormalForm_HH
// ------------------------------------------------------------------------
// $RCSfile: NormalForm.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NormalForm
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include "Algebra/VpsInvMap.h"
#include <complex>

using std::complex;


// Class NormalForm
// ------------------------------------------------------------------------
/// Resonance-free normal form.
//  Find normal form of a truncated Taylor series map. Compute from
//  representation of a map which can be a nil-potent, static or dynamic
//  symplectic map. Implementation of an algorithm described in
//  {CENTER}
//  M. Berz, E. Forest and J. Irwin,{BR}
//  Particle Accelerators, 1989, Vol. 24, pp. 91-107.
//  {/CENTER}

class NormalForm {

public:

    /// Constructor.
    //  Construct normal-form for [b]map[/b].
    NormalForm(const VpsInvMap<double> &map);

    NormalForm();
    NormalForm(const NormalForm &);
    ~NormalForm();

    /// Get number of stable degrees of freedom.
    int degreesOfFreedom() const;

    /// Get normal-form.
    //  Return normal-form map as a Lie transform.
    const Tps<double> &normalForm() const;

    /// Get normalising map.
    //  Return the normalising map as a Lie transform.
    const Tps<double> &normalisingMap() const;

    /// Get eigenvalues.
    //  Return the eigenvalues of the linear part as a complex vector.
    const Vector<complex<double> > &eigenValues() const;

    /// Get eigenvectors.
    //  Return the eigenvectors of the linear part in packed form.
    const Matrix<double> &eigenVectors() const;

    /// Get anharmonicities.
    //  Return the anharmonicities as a symmetric matrix.
    Matrix<double> anharmonicity() const;

    /// Get invariant polynomial.
    //  Return the invariant polynomial for the mode [b]i[/b].
    Tps<double> invariant(int i) const;

protected:

    // Order the modes of the map and associate them to the planes.
    void orderModes(Vector<complex<double> >, Matrix<double>);

private:

    // Not implemented.
    void operator=(const NormalForm &);

    // Representation of the normal form analysis results.
    // ----------------------------------------------------------------------
    // The dimensions of the problem.
    int dimension;

    // Number of degrees of freedom.
    int freedom;

    // The factorised normalising map.
    Tps<double> A_Lie;

    // The factorised normal form map/
    Tps<double> N_Lie;

    // The vector of eigenvalues.
    Vector<complex<double> > lambda;

    // The matrix of eigenvectors.
    Matrix<double> V;
};

#endif // MAD_NormalForm_HH
