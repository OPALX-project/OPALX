#ifndef CLASSIC_Factory_HH
#define CLASSIC_Factory_HH

// ------------------------------------------------------------------------
// $RCSfile: Factory.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Factory
//
// ------------------------------------------------------------------------
// Class category: Construction
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <string>

class ElementBase;
class AttributeSet;


// Class Factory
// ------------------------------------------------------------------------
/// Abstract interface for an element factory.

class Factory {

public:

    Factory();
    virtual ~Factory();

    /// Define a new element.
    //  The element [b]newElement[/b] is linked to the repository.
    //  If an element with the same name exists already, replacement is
    //  rejected, and [b]newElement[/b] is deleted.
    virtual bool define(ElementBase *newElement) = 0;

    /// Erase element by name.
    //  If there is no element with the given [b]name[/b],
    //  the request is ignored.
    virtual void erase(const std::string &name) = 0;

    /// Find element by name.
    //  If an element with the name [b]name[/b] exists,
    //  return a pointer to this element, otherwise return NULL.
    virtual ElementBase *find(const std::string &name) const = 0;

    /// Make new element.
    //  Create a new element with the type [b]type[/b], the name [b]name[/b]
    //  and the attributes in [b]set[/b].  If an element with the name
    //  [b]name[/b] already exists, it is replaced.
    virtual ElementBase *makeElement(const std::string &type, const std::string &name,
                                     const AttributeSet &set) = 0;

    /// Define a new element.
    //  The element [b]newElement[/b] is linked to the repository.
    //  If an element with the same name exists already, it is replaced.
    virtual bool storeElement(ElementBase *newElement) = 0;
};

#endif // CLASSIC_Factory_HH
