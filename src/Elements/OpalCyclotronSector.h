#ifndef OPALX_OPAL_CYCLOTRON_SECTOR_H
#define OPALX_OPAL_CYCLOTRON_SECTOR_H
#include "AbsBeamline/CyclotronSector.h"
#include "AbstractObjects/Definition.h"
#include "Elements/OpalElement.h"

/** @brief Named PSI mirrored trim-coil definition, separate from beamline elements.
 * Supports TYPE="PSI-BFIELD-MIRRORED" only. RMIN/RMAX [m] are shape radii,
 * BMAX [T] is the signed strength, and SLPTC [1/m] is the positive slope.
 * All azimuths are enabled. No legacy millimetre conversion is performed here.
 * @see CyclotronTrimCoil
 */
class OpalTrimCoil : public Definition {
public:
    OpalTrimCoil();
    OpalTrimCoil* clone(const std::string& name) override { return new OpalTrimCoil(name, this); }
    /// Validate the definition without allocating a tracking element.
    void execute() override { model(); }
    bool canReplaceBy(Object* object) override {
        return dynamic_cast<OpalTrimCoil*>(object) != nullptr;
    }
    /** @return Validated parameters by value, suitable for copying into a device view.
     * @throws OpalException For unsupported TYPE, nonfinite values, or invalid radii/slope.
     */
    CyclotronTrimCoil model() const;

private:
    enum { TYPE, RMIN, RMAX, BMAX, SLPTC, SIZE };
    OpalTrimCoil(const std::string& name, OpalTrimCoil* parent) : Definition(name, parent) {}
};

/** @brief Parser and host-side configuration of a placed PSI magnetic sector.
 * FMAPFN identifies the PSI map; SYM defaults to 8. Explicit RMIN/RMAX [m] must
 * match its bounds (1e-10 m comparison tolerance); omitted bounds use the file.
 * VMIN/VMAX default to -0.05/+0.05 m and BSCALE to 1. TRIMCOIL is a string array
 * of named OpalTrimCoil definitions, resolved when update() configures the sector.
 * Common X/Y/Z/THETA/PHI/PSI attributes retain the existing Mode-A placement rules.
 * @see CyclotronSector
 */
class OpalCyclotronSector : public OpalElement {
public:
    OpalCyclotronSector();
    OpalCyclotronSector* clone(const std::string& name) override;
    /** @brief Resolve map/coil definitions and rebuild the occurrence configuration.
     * An empty FMAPFN is tolerated for the parser exemplar; tracking an unconfigured
     * occurrence is rejected by CyclotronSector::initialise().
     * @throws OpalException For invalid symmetry/bounds or an unknown trim-coil name.
     * Map loading and trim-coil validation errors are propagated.
     */
    void update() override;

private:
    enum { FMAPFN = COMMON, SYM, RMIN, RMAX, VMIN, VMAX, BSCALE, TRIMCOIL, SIZE };
    OpalCyclotronSector(const std::string& name, OpalCyclotronSector* parent);
};
#endif
