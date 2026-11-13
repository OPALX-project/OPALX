#include "Elements/OpalCyclotronSector.h"

#include <cmath>
#include <memory>
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
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
CyclotronTrimCoil OpalTrimCoil::model() const {
    CyclotronTrimCoil m{
            Attributes::getReal(itsAttr[RMIN]), Attributes::getReal(itsAttr[RMAX]),
            Attributes::getReal(itsAttr[BMAX]), Attributes::getReal(itsAttr[SLPTC])};
    if (Attributes::getString(itsAttr[TYPE]) != "PSI-BFIELD-MIRRORED" || !(m.rmin > 0)
        || !(m.rmax > m.rmin) || !(m.slope > 0)
        || !std::isfinite(m.rmin + m.rmax + m.slope + m.bmax))
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
    setElement(new CyclotronSector("CYCLOTRONSECTOR"));
}
OpalCyclotronSector::OpalCyclotronSector(const std::string& name, OpalCyclotronSector* parent)
    : OpalElement(name, parent) {
    setElement(new CyclotronSector(name));
}
OpalCyclotronSector* OpalCyclotronSector::clone(const std::string& name) {
    return new OpalCyclotronSector(name, this);
}
std::shared_ptr<ElementBase> OpalCyclotronSector::build(const ElementValues& values) {
    auto element = std::make_shared<CyclotronSector>(values.getName());
    fill(*element, values);
    return element;
}

void OpalCyclotronSector::update() { fill(*getElement(), getValues()); }

void OpalCyclotronSector::fill(ElementBase& element, const ElementValues& values) {
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
        auto* coil = dynamic_cast<OpalTrimCoil*>(OpalData::getInstance()->find(name));
        if (!coil) throw OpalException("OpalCyclotronSector", "Unknown TRIMCOIL: " + name);
        coils.push_back(coil->model());
    }
    static_cast<CyclotronSector*>(&element)->configure(
            map, static_cast<int>(sym), values.getReal("VMIN"), values.getReal("VMAX"),
            values.getReal("BSCALE"), coils);
}
