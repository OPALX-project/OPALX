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
class ElementValues;

namespace opalx::run {

    /// Builds the runtime element of one occurrence from its attribute values.
    using BuildFunction = std::shared_ptr<ElementBase> (*)(const ElementValues& values);

    /**
     * @brief The runtime element type string of an element keyword, e.g. "Multipole" for
     *        QUADRUPOLE.
     * @throws OpalException if @p keyword is not an element keyword.
     */
    std::string elementTypeOf(const std::string& keyword);

    /**
     * @brief Build the runtime element of one line occurrence from its description.
     *
     * The build function of the occurrence's keyword makes the runtime element from the
     * attribute values in @p element. Keywords without a build function yet make an object
     * from the keyword table, copy the attribute values the input set into it, and let its
     * update() fill its runtime element. The parsed input objects are not used.
     *
     * @throws OpalException if the keyword is unknown or an attribute the input set could not
     *         be evaluated or has a type that cannot be copied.
     */
    std::shared_ptr<ElementBase> buildElement(const ElementDescription& element);

    /**
     * @brief Build the runtime elements of every occurrence of @p line, in line order.
     * @note On a RING every element except the #S and #E markers gets the RING membership,
     *       with the RING name as owner, also inside nested LINEs.
     */
    std::vector<std::shared_ptr<ElementBase>> buildElements(const LineDescription& line);

    /**
     * @brief Circumference of a RING: the sum of the arc lengths of its elements, markers left
     *        out.
     */
    double ringLength(const std::vector<std::shared_ptr<ElementBase>>& elements);

}  // namespace opalx::run

#endif  // OPALX_BUILD_ELEMENTS_H
