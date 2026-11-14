#include "Elements/OpalCyclotronSector.h"

#include <cmath>
#include <map>
#include <memory>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Track/FreezeRun.h"
#include "Utilities/OpalException.h"

OpalTrimCoil::OpalTrimCoil()
    : Definition(SIZE, "TRIMCOIL", "Named cyclotron trim-coil field model.") {
    itsAttr[TYPE]  = Attributes::makeString("TYPE", "Trim-coil model.", "PSI-BFIELD-MIRRORED");
    itsAttr[RMIN]  = Attributes::makeReal("RMIN", "Inner shape radius [m].");
    itsAttr[RMAX]  = Attributes::makeReal("RMAX", "Outer shape radius [m].");
    itsAttr[BMAX]  = Attributes::makeReal("BMAX", "Strength [T].");
    itsAttr[SLPTC] = Attributes::makeReal("SLPTC", "Profile slope [1/m].");
    registerOwnership(AttributeHandler::STATEMENT);
}
namespace {
    std::map<std::string, CyclotronTrimCoil>& trimCoils() {
        static std::map<std::string, CyclotronTrimCoil> coils;
        return coils;
    }
}  // namespace
void OpalTrimCoil::execute() { model(ObjectValues(opalx::run::describeObject(*this))); }
void OpalTrimCoil::registerModel(const ObjectValues& values) {
    trimCoils().insert_or_assign(values.getName(), model(values));
}
const CyclotronTrimCoil& OpalTrimCoil::find(const std::string& name) {
    const auto found = trimCoils().find(name);
    if (found == trimCoils().end())
        throw OpalException("OpalCyclotronSector", "Unknown TRIMCOIL: " + name);
    return found->second;
}
CyclotronTrimCoil OpalTrimCoil::model(const ObjectValues& values) {
    CyclotronTrimCoil m{
            values.getReal("RMIN"), values.getReal("RMAX"), values.getReal("BMAX"),
            values.getReal("SLPTC")};
    if (values.getString("TYPE") != "PSI-BFIELD-MIRRORED" || !(m.rmin > 0) || !(m.rmax > m.rmin)
        || !(m.slope > 0) || !std::isfinite(m.rmin + m.rmax + m.slope + m.bmax))
        throw OpalException(
                "OpalTrimCoil",
                "Expected PSI-BFIELD-MIRRORED with finite SI parameters, 0<RMIN<RMAX and SLPTC>0.");
    return m;
}
OpalCyclotronSector::OpalCyclotronSector()
    : OpalElement(SIZE, "CYCLOTRONSECTOR", "PSI median-plane cyclotron magnetic sector.") {
    itsAttr[FMAPFN]   = Attributes::makeString("FMAPFN", "PSI map filename.");
    itsAttr[SYM]      = Attributes::makeReal("SYM", "Number of sectors.", 8);
    itsAttr[RMIN]     = Attributes::makeReal("RMIN", "Map minimum radius [m].");
    itsAttr[RMAX]     = Attributes::makeReal("RMAX", "Map maximum radius [m].");
    itsAttr[VMIN]     = Attributes::makeReal("VMIN", "Minimum vertical coordinate [m].", -0.05);
    itsAttr[VMAX]     = Attributes::makeReal("VMAX", "Maximum vertical coordinate [m].", 0.05);
    itsAttr[BSCALE]   = Attributes::makeReal("BSCALE", "Base-map multiplier.", 1);
    itsAttr[TRIMCOIL] = Attributes::makeStringArray("TRIMCOIL", "Named trim-coil definitions.");
    registerOwnership();
}
std::shared_ptr<ElementBase> OpalCyclotronSector::build(const ObjectValues& values, double /*p0*/) {
    auto element = std::make_shared<CyclotronSector>(values.getName());
    fill(*element, values);
    return element;
}

void OpalCyclotronSector::fill(ElementBase& element, const ObjectValues& values) {
    setCommonFields(element, values);
    const auto filename = values.getString("FMAPFN");
    if (filename.empty()) return;  // no field map: nothing to configure
    auto map         = CyclotronSectorFieldMap::read(filename);
    const double sym = values.getReal("SYM");
    if (!std::isfinite(sym) || sym < 2 || sym > 100000 || sym != std::round(sym))
        throw OpalException("OpalCyclotronSector", "SYM must be an integer >=2.");
    for (const char* a : {"RMIN", "RMAX"}) {
        const double expected =
                std::string(a) == "RMIN" ? map->rmin : map->rmin + (map->nr - 1) * map->dr;
        if (!values.defaultUsed(a)
            && (!std::isfinite(values.getReal(a))
                || std::abs(values.getReal(a) - expected) > 1e-10))
            throw OpalException("OpalCyclotronSector", "RMIN/RMAX must match the field map.");
    }
    std::vector<CyclotronTrimCoil> coils;
    for (const auto& name : values.getStringArray("TRIMCOIL")) {
        coils.push_back(OpalTrimCoil::find(name));
    }
    static_cast<CyclotronSector*>(&element)->configure(
            map, static_cast<int>(sym), values.getReal("VMIN"), values.getReal("VMAX"),
            values.getReal("BSCALE"), coils);
}
