// ------------------------------------------------------------------------
// $RCSfile: ErrorField.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorField
//   Class for the OPAL EFIELD error command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorField.h"
#include "AbstractObjects/BeamSequence.h"
#include "Attributes/Attributes.h"
#include "Errors/Error.h"
#include "Errors/MPHandler.h"
#include "Fields/BMultipoleField.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"
#include <cmath>
#include <iostream>
#include <vector>

using std::cerr;
using std::endl;
using std::max;
using std::vector;


// Class ErrorField
// ------------------------------------------------------------------------

// The attributes of class ErrorField.
namespace {
    enum {
        RADIUS, // The normalising radius in m.
        ORDER,  // The order of component used for normalisation.
        DKR,    // The relative magnitude of error components.
        DK,     // The absolute magnitude of error components in m^(-n).
        ROT,    // The rotation angles in rad.
        SIZE
    };
}


ErrorField::ErrorField():
    Action(SIZE, "EFIELD",
           "The \"EFIELD\" sub-command assigns field errors by magnitude and "
           "rotation angle.") {
    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "Normalising radius in m", 1.0);

    itsAttr[ORDER] = Attributes::makeReal
                     ("ORDER", "Order of component used for normalisation");

    itsAttr[DKR] = Attributes::makeRealArray
                   ("DKR", "Relative magnitude of error components (no dimension)");
    itsAttr[DKR].setDeferred(true);

    itsAttr[DK] = Attributes::makeRealArray
                  ("DK", "Absolute magnitude of error components in m^(-n)");
    itsAttr[DK].setDeferred(true);

    itsAttr[ROT] = Attributes::makeRealArray
                   ("ROT", "Rotation angles in rad");
    itsAttr[ROT].setDeferred(true);

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


ErrorField::ErrorField(const std::string &name, ErrorField *parent):
    Action(name, parent)
{}


ErrorField::~ErrorField()
{}


ErrorField *ErrorField::clone(const std::string &name) {
    return new ErrorField(name, this);
}


void ErrorField::execute() {
    // Warn about redundant itsAttr.
    if(itsAttr[DKR]  &&  itsAttr[DK]) {
        if(Options::warn) {
            cerr << "\n### Warning ### " << *this
                 << "\nBoth relative and absolute errors given; "
                 << "error definitions will use their sum.\n" << endl;
        }
    }

    // Generate the multipole errors.
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

        cerr << endl;
    }
}


void ErrorField::fieldError(const std::string &, int,
                            const BMultipoleField &designField,
                            BMultipoleField &errorField) {
    // Use of relative error requires radius to be given.
    double rad = Attributes::getReal(itsAttr[RADIUS]);
    if(rad == 0.0) rad = 1.0;
    int ord = int(Round(Attributes::getReal(itsAttr[ORDER])));
    double base = sqrt(pow(designField.getNormalComponent(ord + 1), 2) +
                       pow(designField.getSkewComponent(ord + 1), 2));
    if(base == 0.0) base = 1.0;

    // Generate the field error.
    vector<double> relative = Attributes::getRealArray(itsAttr[DKR]);
    vector<double> absolute = Attributes::getRealArray(itsAttr[DK]);
    vector<double> rotation = Attributes::getRealArray(itsAttr[ROT]);
    vector<double>::size_type top =
        max(relative.size(), absolute.size());

    BMultipoleField tempField;
    for(vector<double>::size_type comp = 0; comp < top; comp++) {
        // Compute error amplitude.
        double normal = 0.0;
        double skewed = 0.0;

        if(comp < absolute.size()) {
            normal  = absolute[comp] * absFactor(comp);
        }

        if(comp < relative.size()) {
            normal += relative[comp] * relFactor(comp, ord, rad) * base;
        }

        // Convert to components.
        if(normal != 0.0) {
            if(comp < rotation.size()) {
                double angle = - rotation[comp] * double(comp + 1);
                skewed = normal * sin(angle);
                normal = normal * cos(angle);
            }

            tempField.setNormalComponent(comp + 1, normal);
            tempField.setSkewComponent(comp + 1, skewed);
        }
    }

    if(errorField.order() != 0) {
        oldError++;
        if(Error::block->addError) {
            errorField.addField(tempField);
        } else {
            errorField = tempField;
        }
    } else {
        newError++;
        errorField = tempField;
    }
}