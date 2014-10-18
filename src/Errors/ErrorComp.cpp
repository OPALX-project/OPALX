// ------------------------------------------------------------------------
// $RCSfile: ErrorComp.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorComp
//   Class for OPAL EFCOMP error command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorComp.h"
#include "AbstractObjects/BeamSequence.h"
#include "Attributes/Attributes.h"
#include "Errors/Error.h"
#include "Errors/MPHandler.h"
#include "Fields/BMultipoleField.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"
#include <cmath>
#include <iostream>
#include <vector>

using std::cerr;
using std::max;
using std::vector;


// Class ErrorComp
// ------------------------------------------------------------------------

// The attributes of class ErrorComp.
namespace {
    enum {
        RADIUS,  // The normalising radius in m.
        ORDER,   // The order of component used for normalisation.
        DKNR,    // The relative normal error components.
        DKSR,    // The relative skewed error components.
        DKN,     // The absolute normal error components in m^(-n).
        DKS,     // The absolute skewed error components in m^(-n).
        SIZE
    };
}


ErrorComp::ErrorComp():
    Action(SIZE, "EFCOMP",
           "The \"EFCOMP\" sub-command assigns error fields as multipole "
           "components to elements.") {
    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "Normalising radius in m", 1.0);

    itsAttr[ORDER] = Attributes::makeReal
                     ("ORDER", "Order of component used for normalisation");

    itsAttr[DKNR] = Attributes::makeRealArray
                    ("DKNR", "Relative normal error components (no dimension)");
    itsAttr[DKNR].setDeferred(true);

    itsAttr[DKSR] = Attributes::makeRealArray
                    ("DKSR", "Relative skewed error components (no dimension)");
    itsAttr[DKSR].setDeferred(true);

    itsAttr[DKN] = Attributes::makeRealArray
                   ("DKN", "Absolute normal error components in m^(-n)");
    itsAttr[DKN].setDeferred(true);

    itsAttr[DKS] = Attributes::makeRealArray
                   ("DKS", "Absolute skewed error components in m^(-n)");
    itsAttr[DKS].setDeferred(true);
}


ErrorComp::ErrorComp(const std::string &name, ErrorComp *parent):
    Action(name, parent)
{}


ErrorComp::~ErrorComp()
{}


ErrorComp *ErrorComp::clone(const std::string &name) {
    return new ErrorComp(name, this);
}


void ErrorComp::execute() {
    // Warn about redundant attributes.
    if(itsAttr[DKNR] && itsAttr[DKN]) {
        if(Options::warn) {
            cerr << "\n### Warning ### " << *this
                 << "\nBoth relative and absolute errors given for normal"
                 << " components; error definitions will use their sum.\n"
                 << std::endl;
        }
    }

    if(itsAttr[DKSR] && itsAttr[DKS]) {
        if(Options::warn) {
            cerr << "\n### Warning ### " << *this
                 << "\nBoth relative and absolute errors given for skew"
                 << " components; error definitions will use their sum.\n"
                 << std::endl;
        }
    }

    // Generate the errors.
    oldError = newError = 0;
    MPHandler inserter(*Error::block->itsLine, *this, false);
    inserter.execute();

    if(Options::info) {
        cerr << '\n';

        if(oldError + newError == 0) {
            cerr << "No field errors assigned.\n";
        } else {
            if(oldError != 0) {
                cerr << "Field errors "
                     << (Error::block->addError ? "superposed on " : "replaced on ")
                     << oldError << " element";
                if(oldError > 1) cerr << "s";
                cerr << ".\n";
            }
            if(newError != 0) {
                cerr << "Field errors assigned to " << newError << " element";
                if(newError > 1) cerr << "s";
                cerr << ".\n";
            }
        }

        cerr << std::endl;
    }
}


void ErrorComp::fieldError(const std::string &, int,
                           const BMultipoleField &designField,
                           BMultipoleField &errorField) {
    // Use of relative error requires "RADIUS" to be given.
    double rad = Attributes::getReal(itsAttr[RADIUS]);
    if(rad == 0.0) rad = 1.0;
    int ord = int(Round(Attributes::getReal(itsAttr[ORDER])));
    double base = sqrt(pow(designField.getNormalComponent(ord + 1), 2) +
                       pow(designField.getSkewComponent(ord + 1), 2));
    if(base == 0.0) base = 1.0;

    // Generate the field error.
    vector<double> absNorm = Attributes::getRealArray(itsAttr[DKN]);
    vector<double> absSkew = Attributes::getRealArray(itsAttr[DKS]);
    vector<double> relNorm = Attributes::getRealArray(itsAttr[DKNR]);
    vector<double> relSkew = Attributes::getRealArray(itsAttr[DKSR]);
    vector<double>::size_type top =
        max(max(relNorm.size(), relSkew.size()),
            max(absNorm.size(), absSkew.size()));
    BMultipoleField tempField;

    for(vector<double>::size_type comp = 0; comp < top; comp++) {
        // Factor for relative components.
        double absFac = absFactor(comp);
        double relFac = relFactor(comp, ord, rad) * base;

        double normal = 0.0;
        if(comp < absNorm.size()) {
            normal =  absNorm[comp] * absFac;
        }

        if(comp < relNorm.size()) {
            normal += relNorm[comp] * relFac;
        }

        double skewed = 0.0;
        if(comp < absSkew.size()) {
            skewed  = absSkew[comp] * absFac;
        }

        if(comp < relSkew.size()) {
            skewed += relSkew[comp] * relFac;
        }

        tempField.setNormalComponent(comp + 1, normal);
        tempField.setSkewComponent(comp + 1, skewed);
    }

    // Allow setting of field in constant structure.
    BMultipoleField &eField = const_cast<BMultipoleField &>(errorField);
    if(eField.order() != 0) {
        oldError++;
        if(Error::block->addError) {
            eField.addField(tempField);
        } else {
            eField = tempField;
        }
    } else {
        newError++;
        eField = tempField;
    }
}
