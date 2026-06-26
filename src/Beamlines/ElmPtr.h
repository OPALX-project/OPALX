#ifndef OPALX_ElmPtr_HH
#define OPALX_ElmPtr_HH

// ------------------------------------------------------------------------
// $RCSfile: ElmPtr.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElmPtr
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <memory>
#include <utility>
#include "AbsBeamline/ElementBase.h"

// Class ElmPtr.
// ------------------------------------------------------------------------
/// A section of a beam line.
//  A beam line is built as a list of ElmPtr.

class ElmPtr {
public:
    ElmPtr();
    ElmPtr(const ElmPtr&);
    ElmPtr(ElementBase* elem);
    explicit ElmPtr(std::shared_ptr<ElementBase> elem);
    virtual ~ElmPtr();

    /// Apply visitor.
    //  If any error occurs, this method throws an exception.
    virtual void accept(BeamlineVisitor&) const;

    /// Get the element pointer.
    inline ElementBase* getElement() const;

    /// Set the element pointer.
    inline void setElement(ElementBase*);
    inline void setElement(std::shared_ptr<ElementBase> elem);

protected:
    // The pointer to the element.
    std::shared_ptr<ElementBase> itsElement;
};

inline ElementBase* ElmPtr::getElement() const { return itsElement.get(); }

inline void ElmPtr::setElement(ElementBase* elem) {
    if (elem == itsElement.get()) {
        return;
    }
    if (!elem) {
        itsElement.reset();
        return;
    }
    auto shared = elem->weak_from_this();
    if (!shared.expired()) {
        itsElement = shared.lock();
        return;
    }
    itsElement.reset(elem);
}

inline void ElmPtr::setElement(std::shared_ptr<ElementBase> elem) { itsElement = std::move(elem); }

#endif  // OPALX_ElmPtr_HH
