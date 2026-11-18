
#include "Structure/BinningCmd.h"
#include "Track/FreezeRun.h"

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
    const ObjectValues values(opalx::run::describeObject(*this));
    getParameterType(values);
    getDumpBinsFrequency(values);
    getTablePrintFrequency(values);
}

Inform& BinningCmd::printInfo(Inform& os, const ObjectValues& values) {
    os << "* ************* B I N N I N G ****************************************************** "
       << endl;
    os << "* BINNING        " << values.getName() << '\n'
       << "* MAXBINS        " << getMaxBins(values) << '\n'
       << "* DESIREDWIDTH   " << getDesiredWidth(values) << '\n'
       << "* BINNINGALPHA   " << getBinningAlpha(values) << '\n'
       << "* BINNINGBETA    " << getBinningBeta(values) << '\n'
       << "* ADAPTIVEBINNING " << (getAdaptiveBinning(values) ? "TRUE" : "FALSE") << '\n'
       << "* PARAMETER      " << getParameter(values) << endl;
    if (dumpBinsToFile(values)) {
        os << "* DUMPBINSFILE   " << getDumpBinsFileName(values) << '\n'
           << "* DUMPBINSFREQ   " << getDumpBinsFrequency(values) << endl;
    }
    os << "* TABLEPRINTFREQ " << getTablePrintFrequency(values) << endl;
    os << "* ********************************************************************************** "
       << endl;

    return os;
}

int BinningCmd::getMaxBins(const ObjectValues& values) {
    return static_cast<int>(values.getReal("MAXBINS"));
}

double BinningCmd::getDesiredWidth(const ObjectValues& values) {
    return values.getReal("DESIREDWIDTH");
}

double BinningCmd::getBinningAlpha(const ObjectValues& values) {
    return values.getReal("BINNINGALPHA");
}

double BinningCmd::getBinningBeta(const ObjectValues& values) {
    return values.getReal("BINNINGBETA");
}

bool BinningCmd::getAdaptiveBinning(const ObjectValues& values) {
    return values.getBool("ADAPTIVEBINNING");
}

std::string BinningCmd::getParameter(const ObjectValues& values) {
    return values.getString("PARAMETER");
}

std::string BinningCmd::getDumpBinsFileName(const ObjectValues& values) {
    std::string filename = values.getString("DUMPBINSFILE");
    // The bins are written as JSON; add the extension if the name has none.
    const std::string extension = ".json";
    if (dumpBinsToFile(values)
        && (filename.size() < extension.size()
            || filename.substr(filename.size() - extension.size()) != extension)) {
        filename += extension;
    }
    return filename;
}

int BinningCmd::getDumpBinsFrequency(const ObjectValues& values) {
    int freq = static_cast<int>(values.getReal("DUMPBINSFREQ"));
    if (dumpBinsToFile(values) && freq < 1) {
        throw OpalException(
                "BinningCmd::getDumpBinsFrequency",
                "DUMPBINSFREQ must be >= 1 if dumping bins to a file.");
    }
    return freq;
}

int BinningCmd::getTablePrintFrequency(const ObjectValues& values) {
    const int freq = static_cast<int>(values.getReal("TABLEPRINTFREQ"));
    if (freq < 0) {
        throw OpalException("BinningCmd::getTablePrintFrequency", "TABLEPRINTFREQ must be >= 0.");
    }
    return freq;
}

bool BinningCmd::dumpBinsToFile(const ObjectValues& values) {
    const std::string filename = values.getString("DUMPBINSFILE");
    return !filename.empty() && filename != "NONE";
}

BinningParameter BinningCmd::getParameterType(const ObjectValues& values) {
    static const std::map<std::string, BinningParameter> stringToParam = {
            {"VELOCITYZ", BinningParameter::VELOCITYZ},
            {"POSITIONZ", BinningParameter::POSITIONZ},
            {"PZ", BinningParameter::PZ},
            {"GAMMAZ", BinningParameter::GAMMAZ}};

    const std::string name = getParameter(values);
    auto it                = stringToParam.find(name);
    if (it == stringToParam.end()) {
        throw OpalException(
                "BinningCmd::getParameterType", "Unknown binning PARAMETER \"" + name + "\"");
    }
    return it->second;
}
