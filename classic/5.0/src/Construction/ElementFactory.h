#ifndef CLASSIC_ElementFactory_HH
#define CLASSIC_ElementFactory_HH

// ------------------------------------------------------------------------
// $RCSfile: ElementFactory.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElementFactory
//
// ------------------------------------------------------------------------
// Class category: Construction
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Construction/Factory.h"
#include <map>


// Class ElementFactory
// ------------------------------------------------------------------------
/// Concrete factory class for CLASSIC elements.
//  Defines the beamline element creation via the factory pattern.
//  When the factory is constructed, empty elements are first created and
//  stored.
//  {p}
//  With the makeElement() method, these elements can be cloned and then
//  filled in with from an ElementImage.  The factory also implements an
//  element repository which can store beam lines.


class ElementFactory: public Factory {

public:

    /// Default constructor.
    //  Fills the repository with all standard element definitions.
    ElementFactory();

    virtual ~ElementFactory();

    /// Define a new element.
    //  The element [b]newElement[/b] is linked to the repository.
    //  If an element with the same name exists already, replacement is
    //  rejected, and [b]newElement[/b] is deleted.
    virtual bool define(ElementBase *newElement);

    /// Erase element by name.
    //  If there is no element with the given [b]name[/b],
    //  the request is ignored.
    virtual void erase(const string &name);

    /// Find element by name.
    //  If an element with the name [b]name[/b] exists,
    //  return a pointer to this element, otherwise return NULL.
    virtual ElementBase *find(const string &name) const;

    /// Make new element.
    //  Create a new element with the type [b]type[/b], the name [b]name[/b]
    //  and the attributes in [b]set[/b].
    //  If an element with the name [b]name[/b] already exists, it is replaced.
    virtual ElementBase *makeElement(const string &type,
                                     const string &name,
                                     const AttributeSet &set);

    /// Define a new element.
    //  The element [b]newElement[/b] is linked to the repository.
    //  If an element with the same name exists already, it is replaced.
    virtual bool storeElement(ElementBase *newElement);

private:

    // Not implemented.
    ElementFactory(const ElementFactory &);
    void operator=(const ElementFactory &);

    // The beamline elements are stored in this map.
    typedef std::map<string, ElementBase *, std::less<string> > MapType;
    MapType inventory;
};

#endif // CLASSIC_ElementFactory_HH
