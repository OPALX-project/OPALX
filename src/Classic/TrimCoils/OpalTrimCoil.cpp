#include "TrimCoils/OpalTrimCoil.h"

// #include "Filters/Filters.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"
#include "AbsBeamline/ElementBase.h"

#include "Utility/IpplInfo.h"
#include <cmath>

extern Inform *gmsg;

using namespace Physics;


// Class OpalTrimCoil
// ------------------------------------------------------------------------

// The attributes of class OpalTrimCoil.
namespace {
    enum {
        TYPE,       // The type of trim coil
        COEF,       // 
        BMAX,       //
        RMIN,       //
        RMAX,       //
        SIZE
    };
}

OpalTrimCoil::OpalTrimCoil():
    Definition(SIZE, "TRIMCOIL",
               "The \"TRIMCOIL\" statement defines a trim coil."),
    trimcoil_m(0) {
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Specifies the type of trim coil: PSI-RING");

    itsAttr[COEF] = Attributes::makeRealArray
                    ("COEF", "List of coefficients");

    itsAttr[BMAX] = Attributes::makeReal
                    ("BMAX", "Maximum magnetic field in Tesla.");

    itsAttr[RMIN] = Attributes::makeReal
                    ("RMIN", "Minimum radius in meters.");

    itsAttr[RMAX] = Attributes::makeReal
                    ("RMAX", "Maximum radius in meters.");

    registerOwnership(AttributeHandler::STATEMENT);

    OpalTrimCoil *defTrimCoil = clone("UNNAMED_TRIMCOIL");
    defTrimCoil->builtin = true;

    try {
        defTrimCoil->update();
        OpalData::getInstance()->define(defTrimCoil);
    } catch(...) {
        delete defTrimCoil;
    }
    
    // FIXME
    this->initOpalTrimCoil();
}


OpalTrimCoil::OpalTrimCoil(const std::string &name, OpalTrimCoil *parent):
    Definition(name, parent),
    trimcoil_m(0)
{}


OpalTrimCoil::~OpalTrimCoil() {
    if (trimcoil_m)
        delete trimcoil_m;
}


bool OpalTrimCoil::canReplaceBy(Object *object) {
    // Can replace only by another WAKE.
    return dynamic_cast<OpalTrimCoil *>(object) != 0;
}


OpalTrimCoil *OpalTrimCoil::clone(const std::string &name) {
    return new OpalTrimCoil(name, this);
}


void OpalTrimCoil::execute() {
    update();
}


OpalTrimCoil *OpalTrimCoil::find(const std::string &name) {
    OpalTrimCoil *trimcoil = dynamic_cast<OpalTrimCoil *>(OpalData::getInstance()->find(name));

    if (trimcoil == 0) {
        throw OpalException("OpalTrimCoil::find()", "OpalTrimCoil \"" + name + "\" not found.");
    }
    return trimcoil;
}


void OpalTrimCoil::update() {
    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_TRIMCOIL");
}


void OpalTrimCoil::initOpalTrimCoil() {
    if (trimcoil_m == 0) {
        *gmsg << "* ************* T R I M C O I L ****************************************************" << endl;
        *gmsg << "OpalTrimCoil::initOpalTrimCoilfunction " << endl;
        *gmsg << "* **********************************************************************************" << endl;
        
        std::string type = Util::toUpper(Attributes::getString(itsAttr[TYPE]));
        
        /*
        double bmax = Attributes::getReal(itsAttr[BMAX]);
        double rmin = Attributes::getReal(itsAttr[RMIN]);
        double rmax = Attributes::getReal(itsAttr[RMAX]);
        std::vector<double> coef = Attributes::getRealArray(itsAttr[COEF]);
        */
        
        *gmsg << *this << endl;
    }
}

Inform& OpalTrimCoil::print(Inform &os) const {
    os << "* ************* T R I M C O I L*****************************************************\n"
       << "* TRIMCOIL       " << getOpalName() << '\n'
       << "* TYPE           " << Attributes::getString(itsAttr[TYPE]) << '\n'
       //<< "* COEF           " << Attributes::getReal(itsAttr[NFREQ]) << '\n'
       << "* BMAX           " << Attributes::getReal(itsAttr[BMAX]) << '\n'
       << "* RMIN           " << Attributes::getReal(itsAttr[RMIN]) << '\n'
       << "* RMAX           " << Attributes::getReal(itsAttr[RMAX]) << '\n'
       << "* ********************************************************************************** " << endl;
    return os;
}
