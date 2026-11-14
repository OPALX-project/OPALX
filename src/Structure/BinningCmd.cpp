
#include "Structure/BinningCmd.h"

#include <map>

#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"

BinningCmd::BinningCmd()
    : Definition(
              BINNING::SIZE, "BINNING",
              "The \"BINNING\" statement defines adaptive binning parameters.") {
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

    itsAttr[BINNING::ADAPTIVEBINNING] = Attributes::makeBool(
            "ADAPTIVEBINNING",
            "If true, merge the initial uniform bins using the adaptive binning algorithm. "
            "If false, keep the uniform MAXBINS histogram. Default: true.",
            true);

    itsAttr[BINNING::DUMPBINSFILE] = Attributes::makeString(
            "DUMPBINSFILE",
            "The json file name for dumping bin configuration to a file. Default: \"NONE\" (no "
            "dumping).",
            "NONE");

    itsAttr[BINNING::DUMPBINSFREQ] = Attributes::makeReal(
            "DUMPBINSFREQ", "The frequency in timesteps of dumping bins to a file. Default: 1.", 1);

    itsAttr[BINNING::TABLEPRINTFREQ] = Attributes::makeReal(
            "TABLEPRINTFREQ",
            "The frequency in global timesteps of printing bin statistics to console. "
            "If 0, printing is disabled. Default: 10.",
            10);
}

BinningCmd::~BinningCmd() {}

BinningCmd* BinningCmd::find(const std::string& name) {
    BinningCmd* bc = dynamic_cast<BinningCmd*>(OpalData::getInstance()->find(name));

    if (bc == nullptr) {
        throw OpalException("BinningCmd::find()", "BinningCmd \"" + name + "\" not found.");
    }
    return bc;
}

void BinningCmd::execute() {
    // Check the values when the BINNING statement is read.
    getParameterType();
    getDumpBinsFrequency();
    getTablePrintFrequency();
}

Inform& BinningCmd::printInfo(Inform& os) const {
    os << "* ************* B I N N I N G ****************************************************** "
       << endl;
    os << "* BINNING        " << getOpalName() << '\n'
       << "* MAXBINS        " << getMaxBins() << '\n'
       << "* DESIREDWIDTH   " << getDesiredWidth() << '\n'
       << "* BINNINGALPHA   " << getBinningAlpha() << '\n'
       << "* BINNINGBETA    " << getBinningBeta() << '\n'
       << "* ADAPTIVEBINNING " << (getAdaptiveBinning() ? "TRUE" : "FALSE") << '\n'
       << "* PARAMETER      " << getParameter() << endl;
    if (dumpBinsToFile()) {
        os << "* DUMPBINSFILE   " << getDumpBinsFileName() << '\n'
           << "* DUMPBINSFREQ   " << getDumpBinsFrequency() << endl;
    }
    os << "* TABLEPRINTFREQ " << getTablePrintFrequency() << endl;
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

bool BinningCmd::getAdaptiveBinning() const {
    return Attributes::getBool(itsAttr[BINNING::ADAPTIVEBINNING]);
}

std::string BinningCmd::getParameter() const {
    return Attributes::getString(itsAttr[BINNING::PARAMETER]);
}

std::string BinningCmd::getDumpBinsFileName() const {
    std::string filename = Attributes::getString(itsAttr[BINNING::DUMPBINSFILE]);
    // The bins are written as JSON; add the extension if the name has none.
    const std::string extension = ".json";
    if (dumpBinsToFile()
        && (filename.size() < extension.size()
            || filename.substr(filename.size() - extension.size()) != extension)) {
        filename += extension;
    }
    return filename;
}

int BinningCmd::getDumpBinsFrequency() const {
    int freq = static_cast<int>(Attributes::getReal(itsAttr[BINNING::DUMPBINSFREQ]));
    if (dumpBinsToFile() && freq < 1) {
        throw OpalException(
                "BinningCmd::getDumpBinsFrequency",
                "DUMPBINSFREQ must be >= 1 if dumping bins to a file.");
    }
    return freq;
}

int BinningCmd::getTablePrintFrequency() const {
    const int freq = static_cast<int>(Attributes::getReal(itsAttr[BINNING::TABLEPRINTFREQ]));
    if (freq < 0) {
        throw OpalException("BinningCmd::getTablePrintFrequency", "TABLEPRINTFREQ must be >= 0.");
    }
    return freq;
}

bool BinningCmd::dumpBinsToFile() const {
    const std::string filename = Attributes::getString(itsAttr[BINNING::DUMPBINSFILE]);
    return !filename.empty() && filename != "NONE";
}

BinningParameter BinningCmd::getParameterType() const {
    static const std::map<std::string, BinningParameter> stringToParam = {
            {"VELOCITYZ", BinningParameter::VELOCITYZ},
            {"POSITIONZ", BinningParameter::POSITIONZ},
            {"PZ", BinningParameter::PZ},
            {"GAMMAZ", BinningParameter::GAMMAZ}};

    const std::string name = getParameter();
    auto it                = stringToParam.find(name);
    if (it == stringToParam.end()) {
        throw OpalException(
                "BinningCmd::getParameterType", "Unknown binning PARAMETER \"" + name + "\"");
    }
    return it->second;
}
