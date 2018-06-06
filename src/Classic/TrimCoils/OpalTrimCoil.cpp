#include "TrimCoils/OpalTrimCoil.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "TrimCoils/TrimCoilBFit.h"
#include "TrimCoils/TrimCoilPhaseFit.h"
#include "TrimCoils/TrimCoilMirrored.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"
#include "Utility/IpplInfo.h"

extern Inform *gmsg;


// Class OpalTrimCoil
// ------------------------------------------------------------------------

// The attributes of class OpalTrimCoil.
namespace {
    enum {
        TYPE,       // The type of trim coil
        COEFNUM,    //
        COEFDENOM,  //
        BMAX,       //
        RMIN,       //
        RMAX,       //
        SLPTC,
        SIZE
    };
}

OpalTrimCoil::OpalTrimCoil():
    Definition(SIZE, "TRIMCOIL",
               "The \"TRIMCOIL\" statement defines a trim coil."),
    trimcoil_m(nullptr) {
    itsAttr[TYPE]      = Attributes::makeString
                         ("TYPE", "Specifies the type of trim coil: PSI-BFIELD, PSI-PHASE, PSI-BFIELD-MIRRORED");

    itsAttr[COEFNUM]   = Attributes::makeRealArray
                         ("COEFNUM", "List of polynomial coefficients for the numerator");

    itsAttr[COEFDENOM] = Attributes::makeRealArray
                         ("COEFDENOM", "List of polynomial coefficients for the denominator");

    itsAttr[BMAX]      = Attributes::makeReal
                         ("BMAX", "Maximum magnetic field in Tesla.");

    itsAttr[RMIN]      = Attributes::makeReal
                         ("RMIN", "Minimum radius in millimeters.");

    itsAttr[RMAX]      = Attributes::makeReal
                         ("RMAX", "Maximum radius in millimeters.");

    itsAttr[SLPTC]      = Attributes::makeReal
                         ("SLPTC", "Slopes of the rising edge [1/mm] (for PSI-BFIELD-MIRRORED)");
    

    registerOwnership(AttributeHandler::STATEMENT);

    OpalTrimCoil *defTrimCoil = clone("UNNAMED_TRIMCOIL");
    defTrimCoil->builtin = true;

    try {
        defTrimCoil->update();
        OpalData::getInstance()->define(defTrimCoil);
    } catch(...) {
        delete defTrimCoil;
    }
}


OpalTrimCoil::OpalTrimCoil(const std::string &name, OpalTrimCoil *parent):
    Definition(name, parent),
    trimcoil_m(nullptr)
{}


OpalTrimCoil::~OpalTrimCoil() {
    
}


bool OpalTrimCoil::canReplaceBy(Object *object) {
    // Can replace only by another trim coil.
    return dynamic_cast<OpalTrimCoil *>(object) != nullptr;
}


OpalTrimCoil *OpalTrimCoil::clone(const std::string &name) {
    return new OpalTrimCoil(name, this);
}


void OpalTrimCoil::execute() {
    update();
}


OpalTrimCoil *OpalTrimCoil::find(const std::string &name) {
    OpalTrimCoil *trimcoil = dynamic_cast<OpalTrimCoil *>(OpalData::getInstance()->find(name));

    if (trimcoil == nullptr) {
        throw OpalException("OpalTrimCoil::find()", "OpalTrimCoil \"" + name + "\" not found.");
    }
    return trimcoil;
}


void OpalTrimCoil::update() {
    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_TRIMCOIL");
}


void OpalTrimCoil::initOpalTrimCoil() {
    if (trimcoil_m != nullptr) return;
        
    std::string type = Util::toUpper(Attributes::getString(itsAttr[TYPE]));
        
    double bmax = Attributes::getReal(itsAttr[BMAX]);
    double rmin = Attributes::getReal(itsAttr[RMIN]);
    double rmax = Attributes::getReal(itsAttr[RMAX]);
        
    if (type == "PSI-BFIELD" || type == "PSI-PHASE") {
        std::vector<double> coefnum   = Attributes::getRealArray(itsAttr[COEFNUM]);
        std::vector<double> coefdenom = Attributes::getRealArray(itsAttr[COEFDENOM]);
        if (type == "PSI-BFIELD")
            trimcoil_m = std::unique_ptr<TrimCoilBFit>     (new TrimCoilBFit    (bmax, rmin, rmax, coefnum, coefdenom));
        else // type == "PSI-PHASE"
            trimcoil_m = std::unique_ptr<TrimCoilPhaseFit> (new TrimCoilPhaseFit(bmax, rmin, rmax, coefnum, coefdenom));

    } else if (type == "PSI-BFIELD-MIRRORED") {
        double slope = Attributes::getReal(itsAttr[SLPTC]);
        trimcoil_m = std::unique_ptr<TrimCoilMirrored>     (new TrimCoilMirrored(bmax, rmin, rmax, slope));
    } else {
        WARNMSG(type << " is not a valid trim coil type" << endl);
    }
        
    *gmsg << level3 << *this << endl;
}

Inform& OpalTrimCoil::print(Inform &os) const {
    os << "* ******************************** T R I M C O I L ********************************\n"
       << "* TRIMCOIL       " << getOpalName() << '\n'
       << "* TYPE           " << Attributes::getString(itsAttr[TYPE]) << '\n';
       
    std::string type = Util::toUpper(Attributes::getString(itsAttr[TYPE]));
    if (type == "PSI-BFIELD" || type == "PSI-PHASE") {
        std::vector<double> coefnum = Attributes::getRealArray(itsAttr[COEFNUM]);
        std::stringstream ss;
        for (std::size_t i = 0; i < coefnum.size(); ++i) {
            ss << ((i > 0) ? "+ " : "") << coefnum[i]
               << ((i > 0) ? (" * x^" + std::to_string(i)) : "") << ' ';
        }
        os << "* POLYNOM NUM    " << ss.str() << '\n';
    
        std::vector<double> coefdenom = Attributes::getRealArray(itsAttr[COEFDENOM]);
        ss.str("");
        for (std::size_t i = 0; i < coefdenom.size(); ++i) {
            ss << ((i > 0) ? "+ " : "") << coefdenom[i]
               << ((i > 0) ? (" * x^" + std::to_string(i)) : "") << ' ';
        }
        os << "* POLYNOM DENOM  " << ss.str() << '\n';
    }

    os << "* BMAX           " << Attributes::getReal(itsAttr[BMAX]) << '\n'
       << "* RMIN           " << Attributes::getReal(itsAttr[RMIN]) << '\n'
       << "* RMAX           " << Attributes::getReal(itsAttr[RMAX]) << '\n';
 
    if (Util::toUpper(Attributes::getString(itsAttr[TYPE])) == "PSI-BFIELD-MIRRORED") {
        os << "* SLPTC          " << Attributes::getReal(itsAttr[SLPTC]) << '\n';
    }
    os << "* *********************************************************************************" << endl;
    return os;
}
