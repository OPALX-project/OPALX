//
// Build the runtime elements of a RUN from its run description.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_BUILD_ELEMENTS_H
#define OPALX_BUILD_ELEMENTS_H

#include <memory>
#include <vector>

#include "Track/RunDescription.h"

class ElementBase;

namespace opalx::run {

    /**
     * @brief Build the runtime element of one line occurrence from its description.
     *
     * A new object of the occurrence's keyword is made from the keyword table, the attribute
     * values the input set are copied into it from @p element, and its update() fills its
     * runtime element. The parsed input objects are not used.
     *
     * @throws OpalException if the keyword is unknown or an attribute the input set could not
     *         be evaluated or has a type that cannot be copied.
     */
    std::shared_ptr<ElementBase> buildElement(const ElementDescription& element);

    /**
     * @brief Build the runtime elements of every occurrence of @p line, in line order.
     * @note On a RING every element except the #S and #E markers gets the RING membership,
     *       as Ring::prepareForTracking() gives the parsed occurrences.
     */
    std::vector<std::shared_ptr<ElementBase>> buildElements(const LineDescription& line);

}  // namespace opalx::run

#endif  // OPALX_BUILD_ELEMENTS_H
