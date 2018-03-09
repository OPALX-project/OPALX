#ifndef OPAL_LMDif_HH
#define OPAL_LMDif_HH 1

// ------------------------------------------------------------------------
// $RCSfile: LMDif.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LMDif
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Algebra/Matrix.h"
#include "Match/MatchState.h"

class QRSolver;
template <class T> class Array1D;
template <class T> class Vector;

// Class LMDif
// ------------------------------------------------------------------------
/// The LMDIF command.
//  This class encapsulates a minimisation according to the method used
//  by Levenberg and Marquardt in
//  [center]
//    B. S. Garbow, K. E. Hillstrom, and J. J. More,
//    User Guide for MINPACK-1, ANL 80-74.
//  [/center]
//  Also described in
//  [center]
//    J. F. Bonnans et al., Optimisation Numerique, pp. 70-72.
//    Springer, Berlin, 1997
//  [/center]
//  Algorithm rewritten following the MINPACK package.

class LMDif: public Action {

public:

    /// Exemplar constructor.
    LMDif();

    virtual ~LMDif();

    /// Make clone.
    virtual LMDif *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    LMDif(const LMDif &);
    void operator=(const LMDif &);

    // Clone constructor.
    LMDif(const std::string &name, LMDif *parent);

    // Compute the Jacobian of the set of functions with respect to the
    // set of variables.
    bool findJacobian(Vector<double> &X,           // Current variable values.
                      Vector<double> &F);          // Current function values.

    // Compute the Levenberg-Marquardt parameter p.
    void lmpar(QRSolver &solver,                   // The QR solver object.
               const Array1D<double> &D,           // The diagonal matrix D.
               const Vector<double> &F,            // The function vector.
               Vector<double> &P);                 // The variable step.

    // The step size fraction.
    double fraction;

    // The step size for the LM algorithm.
    double delta;

    // The parameter for the LM algorithm.
    double mu;

    // The Jacobian of the problem.
    Matrix<double> jacobian;
};

#endif // OPAL_LMDif_HH
