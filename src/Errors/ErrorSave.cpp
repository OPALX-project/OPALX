// ------------------------------------------------------------------------
// $RCSfile: ErrorSave.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorSave
//   Class for OPAL ESAVE error command.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/14 07:07:45 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorSave.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbstractObjects/BeamSequence.h"
#include "Attributes/Attributes.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Errors/AlignHandler.h"
#include "Errors/Error.h"
#include "Errors/MPHandler.h"
#include "Fields/BMultipoleField.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using std::endl;
using std::setw;


// Class ErrorSave
// ------------------------------------------------------------------------

// The attributes of class ErrorSav.
namespace {
    enum {
        FNAME,     // The name of file to be written.
        ALIGN,     // If true, print field errors.
        FIELD,     // If true, print alignment errors.
        RADIUS,    // The normalising radius in m.
        ORDER,     // The order of component used for normalisation.
        SIZE
    };
}


ErrorSave::ErrorSave():
    Action(SIZE, "ESAVE",
           "The \"ESAVE\" statement saves the selected imperfections on a file.") {
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to be written", "ESAVE");
    itsAttr[ALIGN] = Attributes::makeBool
                     ("ALIGN", "If true, print alignment errors");
    itsAttr[FIELD] = Attributes::makeBool
                     ("FIELD", "If true, print field errors");
    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "Normalising radius in m", 1.0);
    itsAttr[ORDER] = Attributes::makeReal
                     ("ORDER", "Order of component used for normalisation");
}


ErrorSave::ErrorSave(const string &name, ErrorSave *parent):
    Action(name, parent)
{}


ErrorSave::~ErrorSave()
{}


ErrorSave *ErrorSave::clone(const string &name) {
    return new ErrorSave(name, this);
}


void ErrorSave::execute() {
    // Prepare output stream.
    string file = Attributes::getString(itsAttr[FNAME]);
    os.open(file.c_str());
    if(os.bad()) {
        throw OpalException("ErrorSave::execute()",
                            "Unable to open output stream \"" + file + "\".");
    }

    std::streamsize old_prec = os.precision(12);

    if(Attributes::getBool(itsAttr[ALIGN])) {
        AlignHandler extractor(*Error::block->itsLine, *this, false);
        extractor.execute();
    }

    if(Attributes::getBool(itsAttr[FIELD])) {
        MPHandler extractor(*Error::block->itsLine, *this, false);
        extractor.execute();
    }

    os.precision(old_prec);
}


void ErrorSave::fieldError(const string &name, int occur,
                           const BMultipoleField &designField,
                           BMultipoleField &errorField) {
    if(errorField.order() != 0) {
        double rad = Attributes::getReal(itsAttr[RADIUS]);
        if(rad == 0.0) rad = 1.0;
        int ord = int(Round(Attributes::getReal(itsAttr[ORDER])));
        double base = sqrt(pow(designField.getNormalComponent(ord + 1), 2) +
                           pow(designField.getSkewComponent(ord + 1), 2));
        if(base == 0.0) base = 1.0;
        double factor = relFactor(ord, 0, rad) / base;

	os << "SELECT, CLEAR;\nSELECT, RANGE = "
	   << name << "[" << occur << "], CLASS = " << name << ";\n"
	   << "EFCOMP, RADIUS = " << setw(16) << rad
	   << ", ORDER = " << ord << ",\n  DKNR = {"
	   << setw(16) << factor *errorField.getNormalComponent(1);
	for(int comp = 1; comp < errorField.order(); ++comp) {
	  double factor = relFactor(ord, comp, rad) / base;
	  os << ',';
	  if(comp % 4 == 0) os << "\n          ";
	  os << setw(16) << factor *errorField.getNormalComponent(comp + 1);
	}
	os << "},\n  DKSR = {"
	   << setw(16) << factor *errorField.getSkewComponent(1);
	for(int comp = 1; comp < errorField.order(); ++comp) {
	  double factor = relFactor(ord, comp, rad) / base;
                os << ',';
                if(comp % 4 == 0) os << "\n          ";
                os << setw(16) << factor *errorField.getSkewComponent(comp + 1);
	}
	os << "};" << endl;
    }
}


void ErrorSave::misalignment(const AlignWrapper &wrap, int occur) {
    Euclid3D offset = wrap.offset();
    double dx, dy, ds, vx, vy, vz;

    if(! offset.isIdentity()) {
        const string &name = wrap.getElement()->getName();
        offset.getAll(dx, dy, ds, vx, vy, vz);

	os << "SELECT, CLEAR;\nSELECT, RANGE = "
	   << name << "[" << occur << "], CLASS = " << name
	   << ";\nEALIGN, DELTAX     = " << dx
	   << ", DELTAY     = " << setw(16) << dy
	   << ", DELTAS     = " << setw(16) << ds
	   << ",\n  DTHETA = " << setw(16) << vx
	   << ", DPHI   = " << setw(16) << vy
	   << ", DPSI   = " << setw(16) << vz
	   << ';' << endl;
    }
}
