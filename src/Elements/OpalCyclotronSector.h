#ifndef OPALX_OPAL_CYCLOTRON_SECTOR_H
#define OPALX_OPAL_CYCLOTRON_SECTOR_H
#include "AbsBeamline/CyclotronSector.h"
#include "AbstractObjects/Definition.h"
#include "Elements/OpalElement.h"

/// Named PSI mirrored trim-coil definition; all public attributes use SI units.
class OpalTrimCoil : public Definition {
public:
    OpalTrimCoil();
    OpalTrimCoil* clone(const std::string& name) override { return new OpalTrimCoil(name, this); }
    void execute() override { model(); }
    bool canReplaceBy(Object* object) override {
        return dynamic_cast<OpalTrimCoil*>(object) != nullptr;
    }
    CyclotronTrimCoil model() const;

private:
    enum { TYPE, RMIN, RMAX, BMAX, SLPTC, SIZE };
    OpalTrimCoil(const std::string& name, OpalTrimCoil* parent) : Definition(name, parent) {}
};

/// Parser for a placed PSI cyclotron magnetic sector.
class OpalCyclotronSector : public OpalElement {
public:
    OpalCyclotronSector();
    OpalCyclotronSector* clone(const std::string& name) override;
    void update() override;

private:
    enum { FMAPFN = COMMON, SYM, RMIN, RMAX, VMIN, VMAX, BSCALE, TRIMCOIL, SIZE };
    OpalCyclotronSector(const std::string& name, OpalCyclotronSector* parent);
};
#endif
