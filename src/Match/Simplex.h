#ifndef OPAL_Simplex_HH
#define OPAL_Simplex_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Simplex.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Simplex
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Algebra/Array1D.h"
#include "Match/MatchState.h"

template <class T> class Vector;


// Class Simplex
// ------------------------------------------------------------------------
/// The SIMPLEX command.
//  This class encapsulates a minimisation according to the SIMPLEX method
//  taken from the MINUIT package.

class Simplex: public Action {

public:

    /// Exemplar constructor.
    Simplex();

    virtual ~Simplex();

    /// Make clone.
    virtual Simplex *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Simplex(const Simplex &);
    void operator=(const Simplex &);

    // Clone constructor.
    Simplex(const std::string &name, Simplex *parent);

    // Razzia routine from MINUIT.
    void razzia(double Fnew, Vector<double> &Xnew);

    // The simplex used in the method, and the function values.
    Array1D < Vector<double> > Xsim;
    Array1D <double>           Fsim;

    // The current best and worst points in the simplex.
    int jl, jh;
};

#endif // OPAL_Simplex_HH
