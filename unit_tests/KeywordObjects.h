//
// Make input objects in tests the way the parser does: from the keyword table.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_UNIT_TESTS_KEYWORD_OBJECTS_H
#define OPALX_UNIT_TESTS_KEYWORD_OBJECTS_H

#include <memory>
#include <stdexcept>
#include <string>

#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Track/BuildElements.h"
#include "Track/FreezeRun.h"

namespace opalx::test {

    /**
     * @brief Register keyword @p keyword (class @p T) in the keyword table, if it is not yet.
     * @return The keyword's object, owned by OpalData.
     */
    template <class T>
    T& keyword(const std::string& keyword) {
        OpalData* opal = OpalData::getInstance();
        if (opal->findKeyword(keyword) == nullptr) opal->addKeyword(new T());
        T* model = dynamic_cast<T*>(opal->findKeyword(keyword));
        if (model == nullptr) throw std::logic_error("keyword " + keyword + " has another class");
        return *model;
    }

    /**
     * @brief Make the object an input statement "NAME: KEYWORD;" makes: a new object of
     *        @p keyword with the keyword's attribute defaults.
     * @note The object is not stored in OpalData.
     */
    template <class T>
    std::unique_ptr<T> makeKeywordObject(const std::string& keywordName, const std::string& name) {
        keyword<T>(keywordName);
        std::unique_ptr<Object> object(OpalData::getInstance()->makeObject(keywordName, name));
        if (dynamic_cast<T*>(object.get()) == nullptr)
            throw std::logic_error("keyword " + keywordName + " made no object");
        return std::unique_ptr<T>(static_cast<T*>(object.release()));
    }

    /**
     * @brief An element made from the keyword table and placed at ELEMEDGE = 0.
     *
     * Every input element must be placed (ELEMEDGE or a pose); building refuses one that is
     * not.
     */
    template <class T>
    std::unique_ptr<T> makePlacedElement(const std::string& keywordName, const std::string& name) {
        auto element        = makeKeywordObject<T>(keywordName, name);
        Attribute* elemEdge = element->findAttribute("ELEMEDGE");
        if (elemEdge == nullptr) throw std::logic_error(keywordName + " has no ELEMEDGE");
        Attributes::setReal(*elemEdge, 0.0);
        return element;
    }

    /**
     * @brief Build the runtime element of input object @p object, as a RUN does for an
     *        occurrence of it in the tracked line.
     */
    inline std::shared_ptr<ElementBase> buildFrom(const Object& object) {
        opalx::run::ElementDescription element;
        element.object = opalx::run::describeObject(object);
        return opalx::run::buildElement(element);
    }

}  // namespace opalx::test

#endif  // OPALX_UNIT_TESTS_KEYWORD_OBJECTS_H
