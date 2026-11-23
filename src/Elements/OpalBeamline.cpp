//
// Class OpalBeamline
//   :FIXME: add class description
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Elements/OpalBeamline.h"

#include "AbsBeamline/Monitor.h"
#include "Elements/PlacementResolver.h"

#include "AbstractObjects/OpalData.h"
#include "Moments/BeamMoments.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/MeshGenerator.h"
#include "Utilities/LogicalError.h"
#include "Utilities/Util.h"

#include <filesystem>
#include <fstream>
#include <limits>
#include <regex>

namespace {
    /// Order the element list by ascending field-start position (path length s), ties broken by
    /// name. The field start is the ELEMEDGE for path-length-placed (Mode B) elements and 0 for
    /// 6D-pose (Mode A) ones — getElementPosition() throws when ELEMEDGE is unset, so it is only
    /// read for Mode B. This reproduces the former BeamlineFieldElement::SortAsc.
    double fieldStart(const std::shared_ptr<ElementBase>& e) {
        return e->isElementPositionSet() ? e->getElementPosition() : 0.0;
    }

    bool byFieldStart(
            const std::shared_ptr<ElementBase>& a, const std::shared_ptr<ElementBase>& b) {
        return fieldStart(a) < fieldStart(b)
               || (fieldStart(a) == fieldStart(b) && a->getName() < b->getName());
    }
}  // namespace

OpalBeamline::OpalBeamline()
    : elements_m(), prepared_m(false), referencePathPlacementCompiled_m(false) {}

OpalBeamline::OpalBeamline(const Vector_t<double, 3>& origin, const Quaternion& rotation)
    : elements_m(),
      prepared_m(false),
      referencePathPlacementCompiled_m(false),
      coordTransformationTo_m(origin, rotation) {}

ElementSet OpalBeamline::getElements(const Vector_t<double, 3>& x) {
    ElementSet elementSet;
    const auto selectIfInside = [&](const std::shared_ptr<ElementBase>& element) {
        if (element->isInside(transformToLocalCS(element, x))) elementSet.insert(element);
    };
    // Before placement is final the grid cannot be built; test every element then.
    if (!elementMapping_m && !referencePathPlacementCompiled_m) {
        for (const auto& element : elements_m)
            selectIfInside(element);
        return elementSet;
    }
    for (const auto& element : getElementMapping().elementsAt(x))
        selectIfInside(element);
    return elementSet;
}

void OpalBeamline::buildElementMapping(const std::array<unsigned, 3>& grid, double reach) {
    elementMapping_m = std::make_shared<const GeneralElementMapping>(elements_m, grid, reach);
}

const GeneralElementMapping& OpalBeamline::getElementMapping() const {
    if (!elementMapping_m) {
        throw LogicalError(
                "OpalBeamline::getElementMapping",
                "buildElementMapping() was not called after the last change of the elements.");
    }
    return *elementMapping_m;
}

ElementSet OpalBeamline::getElements() {
    ElementSet elementSet;
    for (auto& item : elements_m) {
        elementSet.insert(item);
    }
    return elementSet;
}

unsigned long OpalBeamline::getFieldAt(
        const unsigned int& /*index*/, const Vector_t<double, 3>& /*pos*/, const long& /*sindex*/,
        const double& /*t*/, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    unsigned long rtv = 0x00;

    return rtv;
}

unsigned long OpalBeamline::getFieldAt(
        const Vector_t<double, 3>& position, const Vector_t<double, 3>& momentum, const double& t,
        Vector_t<double, 3>& Ef, Vector_t<double, 3>& Bf) {
    unsigned long rtv = 0x00;

    ElementSet elements = getElements(position);

    ElementSet::const_iterator it        = elements.begin();
    const ElementSet::const_iterator end = elements.end();

    for (; it != end; ++it) {
        ElementType type = (*it)->getType();
        if (type == ElementType::MARKER || type == ElementType::MONITOR) continue;

        Vector_t<double, 3> localR = transformToLocalCS(*it, position);
        Vector_t<double, 3> localP = rotateToLocalCS(*it, momentum);
        Vector_t<double, 3> localE(0.0), localB(0.0);

        (*it)->applyToReferenceParticle(localR, localP, t, localE, localB);

        Ef += rotateFromLocalCS(*it, localE);
        Bf += rotateFromLocalCS(*it, localB);
    }

    //         if(section.hasWake()) {
    //             rtv |= BEAMLINE_WAKE;
    //         }
    //         if(section.hasParticleMatterInteraction()) {
    //             rtv |= BEAMLINE_PARTICLEMATTERINTERACTION;
    //         }

    return rtv;
}

std::vector<std::array<Vector_t<double, 3>, 2>> OpalBeamline::getReferenceFields(
        const PartBunch_t& bunch) {
    Inform m("OpalBeamline::getReferenceFields");
    const double t = bunch.getT() - 0.5 * bunch.getdT();
    std::vector<std::array<Vector_t<double, 3>, 2>> fields(bunch.getNumParticleContainers());
    for (size_t i = 0; i < fields.size(); ++i) {
        const auto pc = bunch.getParticleContainer(i);
        if (!pc || pc->getTotalNum() == 0) {
            fields[i][0] = Vector_t<double, 3>(0.0);
            fields[i][1] = Vector_t<double, 3>(0.0);
            continue;
        }
        Vector_t<double, 3> externalE(0.0), externalB(0.0);
        getFieldAt(pc->getRefPartR(), pc->getRefPartP(), t, externalE, externalB);
        fields[i][0] = externalB;
        fields[i][1] = externalE;
        m << level5 << "External fields (container " << i << "): externalE = " << externalE
          << ", externalB = " << externalB << endl;
    }
    return fields;
}

void OpalBeamline::switchElementsOff() {
    for (ElementList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++flit)
        (*flit)->goOffline();
}

void OpalBeamline::setMonitorsRecordReference(bool record) {
    for (const auto& element : getElementByType(ElementType::MONITOR))
        static_cast<Monitor*>(element.get())->setRecordReference(record);
}

void OpalBeamline::applyExternalFields(PartBunch_t& bunch) {
    const StepTime time{
            bunch.getT(), bunch.getdT(), static_cast<long long>(bunch.getGlobalTrackStep())};
    forEachElementInBunchFrame(
            bunch,
            [](const ElementBase& element) {
                return element.fieldFrame() == ElementFrame::Local;
            },
            [&time](const std::shared_ptr<ElementBase>& element,
                    const std::shared_ptr<ParticleContainer_t>& pc) {
                element->apply(pc, time);
            });
}

void OpalBeamline::forEachElementInBunchFrame(
        PartBunch_t& bunch, const std::function<bool(const ElementBase&)>& selects,
        const std::function<void(
                const std::shared_ptr<ElementBase>&, const std::shared_ptr<ParticleContainer_t>&)>&
                func) {
    for (size_t ci = 0; ci < bunch.getNumParticleContainers(); ++ci) {
        const auto pc = bunch.getParticleContainer(ci);
        if (!pc) continue;
        const ParticleContainer_t::FrameScope scope(*pc, pc->frame());
        for (const auto& element : candidateElements(bunch, ci)) {
            if (!selects(*element)) continue;
            pc->enterFrame(elementFrame(element, *pc));
            func(element, pc);
        }
    }
}

std::vector<std::shared_ptr<ElementBase>> OpalBeamline::candidateElements(
        const PartBunch_t& bunch, size_t ci) const {
    const auto pc = bunch.getParticleContainer(ci);
    if (!bunch.isPcActive(ci) || !pc) {
        return {};
    }
    // A container without its own or its root's reference particle has no bunch frame.
    const auto root = bunch.getParticleContainer(bunch.isFollower(ci) ? bunch.getPcRoot(ci) : ci);
    if (!root || !root->getReference()) {
        return {};
    }
    pc->requireReferenceFrame("OpalBeamline::candidateElements");

    // Bunch bounds for this container; an empty container queries the point at the reference
    // particle.
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    if (pc->getTotalNum() > 0) {
        const auto bounds = opalx::moments::compute(*pc, opalx::moments::Moment::Bounds);
        rmin              = bounds.minR();
        rmax              = bounds.maxR();
    }

    // The bunch box, extended forward by the distance one step can cover, c dt: a MONITOR
    // records particles that cross its plane during the next step. Elements come in lattice
    // order.
    Vector_t<double, 3> upper = rmax;
    upper(2) += Physics::c * bunch.getdT();
    return getElementMapping().query(rmin, upper, pc->getToLabTrafo());
}

ParticleFrame OpalBeamline::beamlineFrame(const ParticleContainer_t& pc) const {
    return ParticleFrame(ParticleFrame::Kind::Beamline, getCSTrafoLab2Local() * pc.getToLabTrafo());
}

ParticleFrame OpalBeamline::elementFrame(
        const std::shared_ptr<ElementBase>& element, const ParticleContainer_t& pc) const {
    return ParticleFrame(
            ParticleFrame::Kind::Element,
            getMisalignment(element) * (getCSTrafoLab2Local(element) * pc.getToLabTrafo()),
            element.get());
}

size_t OpalBeamline::markBackwardParticlesAtSourcePlane(
        PartBunch_t& bunch, const opalx::spacecharge::DirichletPlaneConfig& plane) const {
    /// \todo this function should probably be integrated as a GunSource element similar to old
    /// OPAL.
    using namespace opalx::spacecharge;
    if (plane.kind == DirichletPlaneType::None) {
        return 0;
    }

    const double sourcePlaneZ = plane.planeZ;
    // Legacy OPAL's SOURCE element is 5 cm long and is shifted upstream from ELEMEDGE.
    // Source::apply deletes only once a particle crosses the element-local entrance plane
    // (Rz <= 0), not when it crosses the cathode/image plane at ELEMEDGE.
    constexpr double legacySourceLength = 0.05;
    const double sourceLossPlaneZ       = sourcePlaneZ - legacySourceLength;

    size_type localTotalMarked = 0;
    const size_t nContainers   = bunch.getNumParticleContainers();
    for (size_t ci = 0; ci < nContainers; ++ci) {
        if (!bunch.isPcActive(ci)) {
            continue;
        }

        auto pc = bunch.getParticleContainer(ci);
        if (!pc || pc->getLocalNum() == 0) {
            continue;
        }

        // Read source-frame coordinates in the kernel; the particles stay where they are.
        const auto toSource = pc->frame().to(beamlineFrame(*pc)).rigid();
        localTotalMarked += pc->markBackwardParticles(toSource, sourceLossPlaneZ);
    }

    size_type globalTotalMarked = 0;
    ippl::Comm->allreduce(localTotalMarked, globalTotalMarked, 1, std::plus<size_type>());
    return static_cast<size_t>(globalTotalMarked);
}

size_t OpalBeamline::applyElementApertures(PartBunch_t& bunch) {
    size_t localMarked = 0;
    for (size_t ci = 0; ci < bunch.getNumParticleContainers(); ++ci) {
        const auto pc = bunch.getParticleContainer(ci);
        for (const auto& element : candidateElements(bunch, ci)) {
            // The check reads local positions in its kernel; the particles stay where they are.
            // Marks are flushed by the single deletion at the end of the step.
            const auto toLocal = pc->frame().to(elementFrame(element, *pc)).rigid();
            localMarked += element->markOutsideAperture(pc, toLocal);
        }
    }

    size_t globalMarked = 0;
    ippl::Comm->allreduce(localMarked, globalMarked, 1, std::plus<size_t>());
    return globalMarked;
}

void OpalBeamline::prepareSections() {
    if (elements_m.empty()) {
        prepared_m = true;
        return;
    }
    elements_m.sort(byFieldStart);
    placeElementsAlongReferencePath();
    elementMapping_m.reset();
    prepared_m = true;
}

void OpalBeamline::print(Inform& /*msg*/) const {}

void OpalBeamline::swap(OpalBeamline& rhs) {
    std::swap(elements_m, rhs.elements_m);
    std::swap(prepared_m, rhs.prepared_m);
    std::swap(referencePathPlacementCompiled_m, rhs.referencePathPlacementCompiled_m);
    std::swap(coordTransformationTo_m, rhs.coordTransformationTo_m);
    std::swap(elementMapping_m, rhs.elementMapping_m);
}

void OpalBeamline::merge(OpalBeamline& rhs) {
    elements_m.insert(elements_m.end(), rhs.elements_m.begin(), rhs.elements_m.end());
    prepared_m                       = false;
    referencePathPlacementCompiled_m = false;
    elementMapping_m.reset();
}

ElementList OpalBeamline::getElementByType(ElementType type) {
    if (type == ElementType::ANY) {
        return elements_m;
    }

    ElementList elements_of_requested_type;
    for (ElementList::iterator fit = elements_m.begin(); fit != elements_m.end(); ++fit) {
        if ((*fit)->getType() == type) {
            elements_of_requested_type.push_back((*fit));
        }
    }
    return elements_of_requested_type;
}

void OpalBeamline::placeElementsAlongReferencePath() {
    if (referencePathPlacementCompiled_m) {
        return;
    }
    // Single PLACE stage: resolve every element's global-to-local transform (6D pose +
    // ELEMEDGE) in one place.
    PlacementResolver::resolve(elements_m, coordTransformationTo_m);
    referencePathPlacementCompiled_m = true;
}

void OpalBeamline::compute3DLattice() { placeElementsAlongReferencePath(); }

void OpalBeamline::save3DLattice() {
    if (ippl::Comm->rank() != 0) return;

    // Write elements in s-sorted order (the order prepareSections established). This is stable
    // for 6D-posed elements too, unlike the former application-order sort.
    elements_m.sort(byFieldStart);

    ElementList::iterator it  = elements_m.begin();
    ElementList::iterator end = elements_m.end();

    std::ofstream pos;
    std::string fileName = Util::combineFilePath(
            {OpalData::getInstance()->getAuxiliaryOutputDirectory(),
             OpalData::getInstance()->getInputBasename() + "_ElementPositions.txt"});
    if (OpalData::getInstance()->getOpenMode() == OpalData::OpenMode::APPEND
        && std::filesystem::exists(fileName)) {
        pos.open(fileName, std::ios_base::app);
    } else {
        pos.open(fileName);
    }

    MeshGenerator mesh;
    for (auto scan = it; scan != end; ++scan) {
        const std::shared_ptr<ElementBase> scanElement = (*scan);
        if (scanElement->getType() == ElementType::DRIFT) {
            continue;
        }

        double minor = 0.0;
        double major = 0.0;
        if (MeshGenerator::getTransverseSupport(*scanElement, minor, major)) {
            mesh.setDriftReference(0.5 * minor, 0.5 * major);
            break;
        }
    }

    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it);
        CoordinateSystemTrafo nominalBody    = getCSTrafoLab2Local(element);
        CoordinateSystemTrafo toBegin        = getNominalEntryTransform(element);
        CoordinateSystemTrafo toEnd          = getNominalExitTransform(element);
        Vector_t<double, 3> entry3D          = toBegin.getOrigin();
        Vector_t<double, 3> exit3D           = toEnd.getOrigin();

        mesh.add(*(element.get()));

        if (element->getGeometry().isBend()) {
            std::vector<Vector_t<double, 3>> designPath = element->getGeometry().getDesignPath();
            unsigned int size                           = designPath.size();

            unsigned int minNumSteps = std::max(
                    20u,
                    static_cast<unsigned int>(std::ceil(
                            std::abs(element->getGeometry().getBendAngle() * Units::rad2deg))));

            unsigned int frequency =
                    std::max(1u, static_cast<unsigned int>(std::floor((double)size / minNumSteps)));

            pos << std::setw(30) << std::left
                << std::string("\"ENTRY EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2) << std::setw(18)
                << std::setprecision(10) << entry3D(0) << std::setw(18) << std::setprecision(10)
                << entry3D(1) << "\n";

            Vector_t<double, 3> position = nominalBody.transformFrom(designPath.front());
            pos << std::setw(30) << std::left
                << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << position(2) << std::setw(18)
                << std::setprecision(10) << position(0) << std::setw(18) << std::setprecision(10)
                << position(1) << std::endl;

            for (unsigned int i = frequency; i + 1 < size; i += frequency) {
                position = nominalBody.transformFrom(designPath[i]);
                pos << std::setw(30) << std::left
                    << std::string("\"MID: ") + element->getName() + std::string("\"")
                    << std::setw(18) << std::setprecision(10) << position(2) << std::setw(18)
                    << std::setprecision(10) << position(0) << std::setw(18)
                    << std::setprecision(10) << position(1) << std::endl;
            }

            position = nominalBody.transformFrom(designPath.back());
            pos << std::setw(30) << std::left
                << std::string("\"END: ") + element->getName() + std::string("\"") << std::setw(18)
                << std::setprecision(10) << position(2) << std::setw(18) << std::setprecision(10)
                << position(0) << std::setw(18) << std::setprecision(10) << position(1)
                << std::endl;

            pos << std::setw(30) << std::left
                << std::string("\"EXIT EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << exit3D(2) << std::setw(18)
                << std::setprecision(10) << exit3D(0) << std::setw(18) << std::setprecision(10)
                << exit3D(1) << std::endl;
        } else {
            pos << std::setw(30) << std::left
                << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2) << std::setw(18)
                << std::setprecision(10) << entry3D(0) << std::setw(18) << std::setprecision(10)
                << entry3D(1) << "\n";

            pos << std::setw(30) << std::left
                << std::string("\"END: ") + element->getName() + std::string("\"") << std::setw(18)
                << std::setprecision(10) << exit3D(2) << std::setw(18) << std::setprecision(10)
                << exit3D(0) << std::setw(18) << std::setprecision(10) << exit3D(1) << std::endl;
        }
    }
    elements_m.sort(byFieldStart);
    mesh.write(OpalData::getInstance()->getInputBasename());
}

namespace {
    std::string parseInput() {
        std::ifstream in(OpalData::getInstance()->getInputFn());
        std::string source("");
        std::string str;
        char testBit;
        const std::string commentFormat("");
        const std::regex empty("^[ \t]*$");
        const std::regex lineEnd(";");
        const std::string lineEndFormat(";\n");
        const std::regex cppCommentExpr("//.*");
        const std::regex cCommentExpr(
                "/\\*.*?\\*/");  // "/\\*(?>[^*/]+|\\*[^/]|/[^*])*(?>(?R)(?>[^*/]+|\\*[^/]|/[^*])*)*\\*/"
        bool priorEmpty = true;

        in.get(testBit);
        // good() also ends the loop if the input file could not be opened.
        while (in.good()) {
            in.putback(testBit);

            std::getline(in, str);
            str = std::regex_replace(
                    str, cppCommentExpr, commentFormat, std::regex_constants::format_default);
            str = std::regex_replace(
                    str, empty, commentFormat, std::regex_constants::format_default);
            if (!str.empty()) {
                source += str;  // + '\n';
                priorEmpty = false;
            } else if (!priorEmpty) {
                source += "##EMPTY_LINE##";
                priorEmpty = true;
            }

            in.get(testBit);
        }

        source = std::regex_replace(source, cCommentExpr, commentFormat);
        source = std::regex_replace(
                source, lineEnd, lineEndFormat,
                std::regex_constants::match_default | std::regex_constants::format_default);

        // Since the positions of the elements are absolute in the laboratory coordinate system we
        // have to make sure that the line command doesn't provide an origin and orientation.
        // Everything after the sequence of elements can be deleted and only "LINE = (...);", the
        // first sub-expression (denoted by '\1'), should be kept.
        const std::regex lineCommand(
                "(LINE[ \t]*=[ \t]*\\([^\\)]*\\))[ \t]*,[^;]*;", std::regex::icase);
        const std::string firstSubExpression("\\1;");
        source = std::regex_replace(source, lineCommand, firstSubExpression);

        return source;
    }

    unsigned int getMinimalSignificantDigits(double num, const unsigned int maxDigits) {
        char buf[32];
        snprintf(buf, 32, "%.*f", maxDigits + 1, num);
        std::string numStr(buf);
        unsigned int length = numStr.length();

        unsigned int numDigits = maxDigits;
        unsigned int i         = 2;
        while (i < maxDigits + 1 && numStr[length - i] == '0') {
            --numDigits;
            ++i;
        }

        return numDigits;
    }

    std::string round2string(double num, const unsigned int maxDigits) {
        char buf[64];

        snprintf(buf, 64, "%.*f", getMinimalSignificantDigits(num, maxDigits), num);

        return std::string(buf);
    }
}  // namespace

void OpalBeamline::save3DInput() {
    if (ippl::Comm->rank() != 0) return;

    ElementList::iterator it  = elements_m.begin();
    ElementList::iterator end = elements_m.end();

    std::string input = parseInput();
    std::string fname = Util::combineFilePath(
            {OpalData::getInstance()->getAuxiliaryOutputDirectory(),
             OpalData::getInstance()->getInputBasename() + "_3D.opal"});
    std::ofstream pos(fname);

    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it);
        std::string elementName              = element->getName();
        const std::regex replacePSI(
                "(" + elementName + "\\s*:[^\\n]*)PSI\\s*=[^,;]*,?", std::regex::icase);
        input = std::regex_replace(input, replacePSI, "\\1\\2");

        const std::regex replaceELEMEDGE(
                "(" + elementName + "\\s*:[^\\n]*)ELEMEDGE\\s*=[^,;]*(.)", std::regex::icase);

        CoordinateSystemTrafo cst  = getCSTrafoLab2Local(element);
        Vector_t<double, 3> origin = cst.getOrigin();
        Vector_t<double, 3> orient =
                Util::getTaitBryantAngles(cst.getRotation().conjugate(), elementName);
        for (unsigned int d = 0; d < 3; ++d)
            orient(d) *= Units::rad2deg;

        std::string x =
                (std::abs(origin(0)) > 1e-10 ? "X = " + round2string(origin(0), 10) + ", " : "");
        std::string y =
                (std::abs(origin(1)) > 1e-10 ? "Y = " + round2string(origin(1), 10) + ", " : "");
        std::string z =
                (std::abs(origin(2)) > 1e-10 ? "Z = " + round2string(origin(2), 10) + ", " : "");

        std::string theta =
                (orient(0) > 1e-10 ? "THETA = " + round2string(orient(0), 6) + " * PI / 180, "
                                   : "");
        std::string phi =
                (orient(1) > 1e-10 ? "PHI = " + round2string(orient(1), 6) + " * PI / 180, " : "");
        std::string psi =
                (orient(2) > 1e-10 ? "PSI = " + round2string(orient(2), 6) + " * PI / 180, " : "");
        std::string coordTrafo = x + y + z + theta + phi + psi;
        if (coordTrafo.length() > 2) {
            coordTrafo = coordTrafo.substr(0, coordTrafo.length() - 2);  // remove last ', '
        }

        std::string position = ("\\1" + coordTrafo + "\\2");

        input = std::regex_replace(input, replaceELEMEDGE, position);
    }

    const std::regex empty("##EMPTY_LINE##");
    const std::string emptyFormat("\n");
    input = std::regex_replace(input, empty, emptyFormat);

    pos << input << std::endl;
}
