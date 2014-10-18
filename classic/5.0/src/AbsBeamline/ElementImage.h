#ifndef CLASSIC_ElementImage_HH
#define CLASSIC_ElementImage_HH

// ------------------------------------------------------------------------
// $RCSfile: ElementImage.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElementImage
//   Contains two std::strings, representing the element name and type, and a
//   map of attributes used to represent this element's state.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/AttributeSet.h"
#include <string>

class ElementBase;


// Class ElementImage
// ------------------------------------------------------------------------
/// An image of an element.
//  Class ElementImage implements an image of an element. It contains two
//  std::strings, the name and the type of the element, and a map of name versus
//  value for all attributes of the element.

class ElementImage: public AttributeSet {

public:

    /// Constructor.
    //  This constructor takes the [b]name[/b] and [b]type[/b] std::strings as
    //  arguments, as well as an AttributeSet mapping attribute names to
    //  values.
    ElementImage(const std::string &name, const std::string &type,
                 const AttributeSet &map);

    ElementImage();
    ElementImage(const ElementImage &);
    virtual ~ElementImage();

    /// Assignment operator.
    const ElementImage &operator=(const ElementImage &);

    /// Set element name.
    void setName(const std::string &name);

    /// Get element name.
    const std::string &getName() const;

    /// Set element type std::string.
    void setType(const std::string &type);

    /// Get element type std::string.
    const std::string &getType() const;

private:

    // String representing element name
    std::string elementName;

    // String representing element type
    std::string elementType;
};

#endif // CLASSIC_ElementImage_HH
