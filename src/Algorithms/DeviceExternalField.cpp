// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include "Algorithms/DeviceExternalField.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"

std::vector<std::shared_ptr<ElementBase>> device_external::orderedCandidates(
        OpalBeamline& beamline, const std::set<std::shared_ptr<ElementBase>>& candidates) {
    auto remaining = candidates;
    std::vector<std::shared_ptr<ElementBase>> ordered;
    ordered.reserve(candidates.size());
    for (const auto& element : beamline.getElementByType(ElementType::ANY)) {
        if (remaining.erase(element)) ordered.push_back(element);
    }
    if (!remaining.empty())
        throw OpalException("DeviceExternalField::orderedCandidates",
                "Candidate element is not a prepared beamline occurrence.");
    return ordered;
}

const char* device_external::Builder::unsupportedReason(const ElementBase& element) {
    // COF uses nominal placement. Do not select a controller which would omit
    // a particle-only misalignment from its support queries.
    const auto misalignment = element.getMisalignment();
    const auto origin = misalignment.getOrigin();
    const auto rotation = misalignment.getRotationMatrix();
    for (unsigned i = 0; i < 3; ++i) {
        if (origin(i) != 0) return "Misalignments are not supported by device boundary control: ";
        for (unsigned j = 0; j < 3; ++j)
            if (rotation(i, j) != (i == j ? 1. : 0.))
                return "Misalignments are not supported by device boundary control: ";
    }
    if (const auto* cavity = dynamic_cast<const VariableRFCavity*>(&element)) {
        return std::isfinite(cavity->getWidth()) && cavity->getWidth() > 0
                       && std::isfinite(cavity->getHeight()) && cavity->getHeight() > 0
                       && std::isfinite(cavity->getLength()) && cavity->getLength() > 0
                ? nullptr : "Device RF boundary control requires finite positive dimensions: ";
    }
    const auto higherBendComponents = [](const auto& magnet) {
        for (int n = 2; n < magnet.maxNormal_m; ++n)
            if (magnet.normalComponentsHost_m(n) != 0) return true;
        for (int n = 2; n < magnet.maxSkew_m; ++n)
            if (magnet.skewComponentsHost_m(n) != 0) return true;
        return false;
    };
    bool higher = false;
    switch (element.getType()) {
        case ElementType::DRIFT:
        case ElementType::MARKER:
        case ElementType::MONITOR:
            return nullptr;
        case ElementType::MULTIPOLE: {
            const auto& magnet = dynamic_cast<const Multipole&>(element);
            for (size_t n = 2; n < magnet.getMaxNormalComponentIndex(); ++n)
                higher = higher || magnet.getNormalComponent(n) != 0;
            for (size_t n = 2; n < magnet.getMaxSkewComponentIndex(); ++n)
                higher = higher || magnet.getSkewComponent(n) != 0;
            break;
        }
        case ElementType::SBEND:
            higher = higherBendComponents(dynamic_cast<const SBend&>(element));
            break;
        case ElementType::RBEND:
            higher = higherBendComponents(dynamic_cast<const RBend&>(element));
            break;
        default:
            return "Device boundary control has no descriptor for ";
    }
    return higher ? "Higher multipoles are not supported by device boundary control: " : nullptr;
}

bool device_external::Builder::supports(OpalBeamline& beamline) {
    for (const auto& element : beamline.getElements())
        if (unsupportedReason(*element)) return false;
    return true;
}

device_external::Lattice device_external::Builder::build(OpalBeamline& beamline) {
    const auto source = beamline.getElements();
    Kokkos::View<Element*> data("Ring::fieldDescriptors", source.size());
    auto host = Kokkos::create_mirror_view(data);
    Lattice result;
    size_t index = 0;
    for (const auto& element : source) {
        auto& out = host(index++);
        if (const auto* reason = unsupportedReason(*element))
            throw OpalException("DeviceExternalField", std::string(reason) + element->getName());
        const auto frame = beamline.getCSTrafoLab2Local(element);
        out.frame.origin = frame.getOrigin();
        out.frame.rotation = frame.getRotationMatrix();
        element->getFieldExtent(out.begin, out.end);
        const auto aperture = element->getAperture();
        out.aperture = aperture.first;
        out.apertureX = aperture.second.at(0);
        out.apertureY = aperture.second.at(1);
        if (const auto* cavity = dynamic_cast<const VariableRFCavity*>(element.get())) {
            out.kind = Element::UniformRF;
            out.apertureX = 0.5 * cavity->getWidth();
            out.apertureY = 0.5 * cavity->getHeight();
            result.magneticOnly = false;
        } else switch (element->getType()) {
            case ElementType::DRIFT:
            case ElementType::MARKER:
            case ElementType::MONITOR:
                break;
            case ElementType::MULTIPOLE: {
                const auto& magnet = dynamic_cast<const Multipole&>(*element);
                out.kind = Element::Multipole;
                out.bend.dipoleNormal = magnet.getNormalComponent(0);
                out.bend.quadNormal = magnet.getNormalComponent(1);
                out.bend.dipoleSkew = magnet.getSkewComponent(0);
                out.bend.quadSkew = magnet.getSkewComponent(1);
                break;
            }
            case ElementType::SBEND:
                out.kind = Element::SectorBend;
                out.bend = dynamic_cast<const SBend&>(*element).makeFieldInputs();
                break;
            case ElementType::RBEND:
                out.kind = Element::RectangularBend;
                out.bend = dynamic_cast<const RBend&>(*element).makeFieldInputs();
                break;
            default:
                throw OpalException("DeviceExternalField", "Device boundary control has no descriptor for " + element->getName());
        }
        if (element->getType() != ElementType::MARKER && element->getType() != ElementType::MONITOR) {
            for (double length : {std::abs(out.end - out.begin), element->getGeometry().getArcLength()})
                if (length > 0) result.maximumStep = std::min(result.maximumStep, length / (4 * Physics::c));
        }
    }
    Kokkos::deep_copy(data, host);
    result.elements = data;
    return result;
}
