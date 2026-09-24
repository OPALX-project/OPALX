#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"

namespace endfieldmodel {

AsymmetricEnge::AsymmetricEnge(
        const std::vector<double> aStart, double x0Start, double lambdaStart,
        const std::vector<double> aEnd, double x0End, double lambdaEnd) {
    config_m.engeStart_m.a_m = Enge::makeView(aStart, "AsymmetricEngeStart");
    config_m.engeStart_m.x0_m = x0Start;
    config_m.engeStart_m.lambda_m = lambdaStart;
    // x0 is held in this
    config_m.engeEnd_m.a_m = Enge::makeView(aEnd, "AsymmetricEngeEnd");
    config_m.engeEnd_m.x0_m = x0End;
    config_m.engeEnd_m.lambda_m = lambdaEnd;
    setMaximumDerivative(EngeConfig::max_derivative);
}

void AsymmetricEnge::rescale(double scaleFactor) {
    config_m.engeStart_m.x0_m *= scaleFactor;
    config_m.engeStart_m.lambda_m *= scaleFactor;
    config_m.engeEnd_m.x0_m *= scaleFactor;
    config_m.engeEnd_m.lambda_m *= scaleFactor;
}

std::ostream& AsymmetricEnge::print(std::ostream& out) const {
    out << "AsymmetricEnge start ";
    auto conf = config_m.engeStart_m;
    out << "AsymmetricEnge start function l=" << conf.lambda_m << " x0="
        << conf.x0_m << " c=";
    for (size_t ai = 0; ai < conf.a_m.extent(0); ++ai) {
        out << conf.a_m(ai) << " ";
    }
    conf = config_m.engeStart_m;
    out << "               end   function l=" << conf.lambda_m << " x0="
        << conf.x0_m << " c=";
    for (size_t ai = 0; ai < conf.a_m.extent(0); ++ai) {
        out << conf.a_m(ai) << " ";
    }
    return out;
}
}  // namespace endfieldmodel
