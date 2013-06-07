#ifndef OPAL_Micado_HH
#define OPAL_Micado_HH

// ------------------------------------------------------------------------
// $RCSfile: Micado.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Micado
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:45 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/CorrectionBase.h"
#include "Algorithms/AbstractMapper.h"
#include "MemoryManagement/OwnPtr.h"

template <class T> class Matrix;
template <class T> class Vector;


/// Class Micado
// ------------------------------------------------------------------------
/// OPAL command MICADO.

class Micado: public CorrectionBase {

public:

    /// Exemplar constructor.
    Micado();

    virtual ~Micado();

    /// Make clone.
    virtual Micado *clone(const string &name);

    /// Check validity of the table definition.
    virtual void execute();

private:

    /// The additional attributes for the "Micado" command.
    enum {
        METHOD = CorrectionBase::SIZE, // The method to be used.
        TOL,         // The tolerance for the closed orbit.
        ITERATIONS,  // The number of iterations.
        CORRECTORS,  // The number of correctors to be used.
        PLANE,       // The plane(s) to be treated.
        LISTC1,      // List the correctors before correction.
        LISTC,       // List the correctors during correction.
        LISTC2,      // List the correctors after correction.
        LISTM1,      // List the monitors before correction.
        LISTM2,      // List the monitors after correction.
        SIZE
    };

private:

    // Not implemented.
    Micado(const Micado &);
    void operator=(const Micado &);

    /// Clone constructor.
    Micado(const string &name, Micado *parent);

    // Set the correctors for this plane.
    void applyCorrections(int mode, Vector<double> &X);

    // Find the closed orbit.
    void findClosedOrbit();

    // Set up the influence matrix for one plane.
    void setupInfluence(int mode, Matrix<double> &);

    // Set up the vector of readings for one plane.
    void setupReadings(int mode, Vector<double> &B);

    // Solve the system of equations.
    void solve(int plane, Matrix<double> &A, Vector<double> &B);


    // The mapper being used for the orbit.
    OwnPtr<AbstractMapper> itsMapper;
};

#endif // OPAL_Micado_HH
