#ifndef OPALX_ABS_BEAMLINE_ELEMENT_TYPE_TABLE_H
#define OPALX_ABS_BEAMLINE_ELEMENT_TYPE_TABLE_H

#include "AbsBeamline/ElementBase.h"

namespace opalx {
namespace viz {

    /**
     * @brief One entry of the element-type taxonomy shared by the in-situ
     *        visualization pipeline (JSON element descriptions and field
     *        sampling).
     */
    struct ElementTypeEntry {
        ElementType type;      ///< Enumerator from ElementBase.h.
        const char* name;      ///< ALL-CAPS machine name (JSON + Python).
        bool carriesField;     ///< Element produces an external field.
    };

    /// The single shared taxonomy: element type -> machine name + field flag.
    constexpr ElementTypeEntry elementTypeTable[] = {
        {ElementType::ANY, "ANY", false},
        {ElementType::BEAMLINE, "BEAMLINE", true},
        {ElementType::COLLIMATOR, "COLLIMATOR", false},
        {ElementType::DRIFT, "DRIFT", false},
        {ElementType::LASER, "LASER", false},
        {ElementType::MARKER, "MARKER", false},
        {ElementType::MONITOR, "MONITOR", false},
        {ElementType::PROBE, "PROBE", false},
        {ElementType::SOURCE, "SOURCE", false},
        {ElementType::VACUUM, "VACUUM", false},
        {ElementType::MULTIPOLE, "MULTIPOLE", true},
        {ElementType::MULTIPOLET, "MULTIPOLET", true},
        {ElementType::RFCAVITY, "RFCAVITY", true},
        {ElementType::TRAVELINGWAVE, "TRAVELINGWAVE", true},
        {ElementType::SBEND, "SBEND", true},
        {ElementType::RBEND, "RBEND", true},
        {ElementType::RBEND3D, "RBEND3D", true},
        {ElementType::SOLENOID, "SOLENOID", true},
        {ElementType::CONSTANTEFIELDCAVITY, "CONSTANTEFIELDCAVITY", true},
        {ElementType::CONSTANTFOCUSING, "CONSTANTFOCUSING", true}};

    /// @return ALL-CAPS machine name for a type, or "OTHER" when absent.
    inline const char* machineName(ElementType type) {
        for (const auto& entry : elementTypeTable)
            if (entry.type == type)
                return entry.name;
        return "OTHER";
    }

    /// @return True when the element type produces an external field;
    ///         unknown types default to true (conservative).
    inline bool carriesField(ElementType type) {
        for (const auto& entry : elementTypeTable)
            if (entry.type == type)
                return entry.carriesField;
        return true;
    }

}  // namespace viz
}  // namespace opalx

#endif  // OPALX_ABS_BEAMLINE_ELEMENT_TYPE_TABLE_H
