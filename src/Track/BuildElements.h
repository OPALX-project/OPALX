//
// Build the runtime elements of a RUN from its run description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_BUILD_ELEMENTS_H
#define OPALX_BUILD_ELEMENTS_H

#include <memory>
#include <string>
#include <vector>

#include "Track/RunDescription.h"

class ElementBase;
class ObjectValues;

namespace opalx::run {

    /// Builds the runtime element of one occurrence from its attribute values. @p p0 is the
    /// input variable P0 (eV/c), which scales the magnet strengths K0, K1, ... and KS.
    using BuildFunction = std::shared_ptr<ElementBase> (*)(const ObjectValues& values, double p0);

    /**
     * @brief The runtime element type string of an element keyword, e.g. "Multipole" for
     *        QUADRUPOLE.
     * @throws OpalException if @p keyword is not an element keyword.
     */
    std::string elementTypeOf(const std::string& keyword);

    /**
     * @brief Build the runtime element of one line occurrence from its description.
     *
     * The build function of the occurrence's keyword (OpalX::build()) makes the runtime
     * element from the attribute values in @p element. The parsed input objects are not used.
     *
     * @param p0 The input variable P0 (eV/c) of the run (RunContext::p0).
     * @throws OpalException if the keyword is unknown or an attribute the input set could not
     *         be evaluated, or the attribute values are not valid for the element.
     */
    std::shared_ptr<ElementBase> buildElement(const ElementDescription& element, double p0);

    /**
     * @brief Build the runtime elements of every occurrence of @p line, in line order.
     * @note On a RING every element except the #S and #E markers gets the RING membership,
     *       with the RING name as owner, also inside nested LINEs.
     */
    std::vector<std::shared_ptr<ElementBase>> buildElements(const LineDescription& line, double p0);

    /// True if @p keyword defines a time dependence, an end-field model (ENGE, TANH, ...) or a
    /// cyclotron trim coil (TRIMCOIL).
    bool isModelKeyword(const std::string& keyword);

    /**
     * @brief Make each time dependence, end-field model and trim coil of @p models and register
     *        it under its name (AbstractTimeDependence, EndFieldModelManager, OpalTrimCoil),
     *        where the elements look it up. buildElements() calls it for the models of the line.
     */
    void registerModels(const std::vector<ObjectDescription>& models);

    /**
     * @brief Circumference of a RING: the sum of the arc lengths of its elements, markers left
     *        out.
     */
    double ringLength(const std::vector<std::shared_ptr<ElementBase>>& elements);

}  // namespace opalx::run

#endif  // OPALX_BUILD_ELEMENTS_H
