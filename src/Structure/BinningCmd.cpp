
#include "Structure/BinningCmd.h"

#include <map>

#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"

BinningCmd::BinningCmd()
    : Definition(
          BINNING::SIZE, "BINNING",
          "The \"BINNING\" statement defines adaptive binning parameters."),
      parameterName_m("VELOCITYZ"),
      parameterType_m(BinningParameter::VELOCITYZ) {
    itsAttr[BINNING::MAXBINS] = Attributes::makeReal(
        "MAXBINS",
        "The maximum number of bins used for adaptive binning. "
        "Default: 128.",
        128);

    itsAttr[BINNING::DESIREDWIDTH] = Attributes::makeReal(
        "DESIREDWIDTH",
        "A bias [0, 1] that tries to steer the bin size to the given variable. "
        "Default: 0.1.",
        0.1);

    itsAttr[BINNING::BINNINGALPHA] = Attributes::makeReal(
        "BINNINGALPHA",
        "A value [0, 1] that determines how aggressive the algorithm tries to "
        "reduce the number of bins. Default: 1.0.",
        1.0);

    itsAttr[BINNING::BINNINGBETA] = Attributes::makeReal(
        "BINNINGBETA",
        "A value [0, 1] that determines how aggressive the algorithm tries to "
        "use wider bins. Default: 1.5.",
        1.5);

    itsAttr[BINNING::PARAMETER] = Attributes::makePredefinedString(
        "PARAMETER", "The bunch attribute used for binning.",
        {"VELOCITYZ", "POSITIONZ", "PZ", "GAMMAZ"}, "VELOCITYZ");
}

BinningCmd::BinningCmd(const std::string& name, BinningCmd* parent)
    : Definition(name, parent),
      parameterName_m(parent->parameterName_m),
      parameterType_m(parent->parameterType_m) {}

BinningCmd::~BinningCmd() {}

BinningCmd* BinningCmd::clone(const std::string& name) { return new BinningCmd(name, this); }

BinningCmd* BinningCmd::find(const std::string& name) {
    BinningCmd* bc = dynamic_cast<BinningCmd*>(OpalData::getInstance()->find(name));

    if (bc == nullptr) {
        throw OpalException("BinningCmd::find()", "BinningCmd \"" + name + "\" not found.");
    }
    return bc;
}

void BinningCmd::execute() {
    setParameterType();
    update();
}

void BinningCmd::update() {
    if (itsAttr[BINNING::PARAMETER]) {
        parameterName_m = Attributes::getString(itsAttr[BINNING::PARAMETER]);
    } else {
        parameterName_m = "VELOCITYZ";
    }
}

void BinningCmd::setParameterType() {
    // Ensure we work from the current attribute value.
    parameterName_m = getParameter();

    static const std::map<std::string, BinningParameter> stringToParam = {
        {"VELOCITYZ", BinningParameter::VELOCITYZ},
        {"POSITIONZ", BinningParameter::POSITIONZ},
        {"PZ", BinningParameter::PZ},
        {"GAMMAZ", BinningParameter::GAMMAZ}};

    auto it = stringToParam.find(parameterName_m);
    if (it == stringToParam.end()) {
        throw OpalException(
            "BinningCmd::setParameterType",
            "Unknown binning PARAMETER \"" + parameterName_m + "\"");
    }

    parameterType_m = it->second;
}

Inform& BinningCmd::printInfo(Inform& os) const {
    os << "* ************* B I N N I N G ****************************************************** "
       << endl;
    os << "* BINNING      " << getOpalName() << '\n'
       << "* MAXBINS      " << getMaxBins() << '\n'
       << "* DESIREDWIDTH " << getDesiredWidth() << '\n'
       << "* BINNINGALPHA " << getBinningAlpha() << '\n'
       << "* BINNINGBETA  " << getBinningBeta() << '\n'
       << "* PARAMETER    " << parameterName_m << endl;
    os << "* ********************************************************************************** "
       << endl;

    return os;
}

int BinningCmd::getMaxBins() const {
    return static_cast<int>(Attributes::getReal(itsAttr[BINNING::MAXBINS]));
}

double BinningCmd::getDesiredWidth() const {
    return Attributes::getReal(itsAttr[BINNING::DESIREDWIDTH]);
}

double BinningCmd::getBinningAlpha() const {
    return Attributes::getReal(itsAttr[BINNING::BINNINGALPHA]);
}

double BinningCmd::getBinningBeta() const {
    return Attributes::getReal(itsAttr[BINNING::BINNINGBETA]);
}

std::string BinningCmd::getParameter() {
    return Attributes::getString(itsAttr[BINNING::PARAMETER]);
}

BinningParameter BinningCmd::getParameterType() const { return parameterType_m; }
