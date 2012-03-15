#ifndef OPAL_Period_HH
#define OPAL_Period_HH
// ------------------------------------------------------------------------
// $RCSfile: Period.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Period
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/Twiss.h"
#include "FixedAlgebra/FVector.h"


// Class Period
// ------------------------------------------------------------------------
/// The TWISS command.

class Period: public Twiss {

public:

    /// Exemplar constructor.
    Period();

    virtual ~Period();

    /// Make clone.
    virtual Period *clone(const string &name);

    /// Fill the buffer using the defined algorithm.
    virtual void fill();

    /// Write TFS file for this table.
    virtual void makeTFS(std::ostream &, const CellArray &) const;

    /// Print the table on an ASCII stream.
    virtual void printTable(std::ostream &, const CellArray &) const;

    /// Write TFS descriptors for a TWISS table.
    //  Writes the line and beam info, as well as length and tunes.
    virtual void tfsTwissDescriptors(std::ostream &) const;

private:

    // Not implemented.
    Period(const Period &);
    void operator=(const Period &);

    // Clone constructor.
    Period(const string &name, Period *parent);

    // Find the closed orbit.
    void findClosedOrbit();

    // The initial closed orbit.
    FVector<double, 6> fixPoint;

    // Additional values for class Period.
    enum {
        MICADO = Twiss::SIZE,  // Number of iterations for MICADO
        CORRECTORS,            // Number of correctors for MICADO
        THREAD,                // Name of threader method
        TOLQ,                  // Tolerances for closed orbit search
        TOLP,

        // Computed values (read-only):
        CIRCUM,                // Machine circumference or line length
        Q1,                    // Tunes
        Q2,
        Q3,
        FREQ,                  // Revolution frequencey in Hz.
        U0,                    // Energy loss per turn in MeV
        J1,                    // Damping partition numbers
        J2,
        J3,
        DELTAP,                // Differential momentum variation
        SIZE
    };
};

#endif // OPAL_Period_HH
