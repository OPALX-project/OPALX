// ------------------------------------------------------------------------
// $RCSfile: ErrorPrint.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorPrint
//   Class for OPAL EPRINT error command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorPrint.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbstractObjects/BeamSequence.h"
#include "Attributes/Attributes.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Errors/AlignHandler.h"
#include "Errors/Error.h"
#include "Errors/MPHandler.h"
#include "Fields/BMultipoleField.h"
#include "Utilities/OpalException.h"
#include "Utilities/Round.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using std::endl;
using std::setw;


// Class ErrorPrint
// ------------------------------------------------------------------------

// The attributes of class ErrorPrint.
namespace {
    enum {
        FNAME,     // The name of file to be written.
        ALIGN,     // If true, print field errors.
        FIELD,     // If true, print alignment errors.
        RADIUS,
        ORDER,
        SIZE
    };
}


ErrorPrint::ErrorPrint():
    Action(SIZE, "EPRINT",
           "The \"EPRINT\" statement prints the selected errors on a file.") {
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to be written", "EPRINT");
    itsAttr[ALIGN] = Attributes::makeBool
                     ("ALIGN", "If true, print alignment errors");
    itsAttr[FIELD] = Attributes::makeBool
                     ("FIELD", "If true, print field errors");
    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "Normalising radius in m", 1.0);
    itsAttr[ORDER] = Attributes::makeReal
                     ("ORDER", "Order of component used for normalisation");

    registerOwnership(AttributeHandler::STATEMENT);
}


ErrorPrint::ErrorPrint(const std::string &name, ErrorPrint *parent):
    Action(name, parent)
{}


ErrorPrint::~ErrorPrint()
{}


ErrorPrint *ErrorPrint::clone(const std::string &name) {
    return new ErrorPrint(name, this);
}


void ErrorPrint::execute() {
    // Prepare output stream.
    std::string file = Attributes::getString(itsAttr[FNAME]);
    os.open(file.c_str());
    if(os.bad()) {
        throw OpalException("ErrorPrint::execute()",
                            "Unable to open output stream \"" + file + "\".");
    }

    std::streamsize old_prec = os.precision(12);

    if(Attributes::getBool(itsAttr[ALIGN])) {
        os << "Misalignment errors for selected elements:" << endl;
        AlignHandler extractor(*Error::block->itsLine, *this, false);
        extractor.execute();
    }

    if(Attributes::getBool(itsAttr[FIELD])) {
        os << "Multipole errors for selected elements:" << endl;
        MPHandler extractor(*Error::block->itsLine, *this, false);
        extractor.execute();
    }

    os.precision(old_prec);
}


void ErrorPrint::fieldError(const std::string &name, int occur,
                            const BMultipoleField &designField,
                            BMultipoleField &errorField) {
    os << "Element name \"" << name << "\" , occurrence = " << occur
       << ':' << endl;
    double rad = Attributes::getReal(itsAttr[RADIUS]);
    if(rad == 0.0) rad = 1.0;
    int ord = int(Round(Attributes::getReal(itsAttr[ORDER])));
    double base = sqrt(pow(designField.getNormalComponent(ord + 1), 2) +
                       pow(designField.getSkewComponent(ord + 1), 2));
    if(base == 0.0) base = 1.0;

    for(int comp = 0; comp < errorField.order(); ++comp) {
        double factor = relFactor(ord, comp, rad) / base;
        os << setw(5) << comp << "  "
           << setw(16) << factor *errorField.getNormalComponent(comp + 1)
           << "  " << setw(16)
           << factor *errorField.getSkewComponent(comp + 1) << endl;
    }
}


void ErrorPrint::misalignment(const AlignWrapper &wrap, int occur) {
    Euclid3D offset = wrap.offset();

    const std::string &name = wrap.getElement()->getName();
    double dx, dy, dz, vx, vy, vz;
    offset.getAll(dx, dy, dz, vx, vy, vz);
    os << "Element name \"" << name << "\", occur = " << occur << '\n'
       << setw(12) << dx << ' ' << setw(12) << dy << ' '
       << setw(12) << dz << ' ' << setw(12) << vx << ' '
       << setw(12) << vy << ' ' << setw(12) << vz << endl;
}