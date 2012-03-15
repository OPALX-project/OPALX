#ifndef OPAL_OpalElement_HH
#define OPAL_OpalElement_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalElement.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalElement
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:32:23 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Element.h"
#include "Elements/AttCell.h"
#include "MemoryManagement/OwnPtr.h"
#include <map>

class Statement;


// Class OpalElement
// ------------------------------------------------------------------------
/// Base class for all beam line elements.
//  This class factors out all special behaviour for the DOOM interface
//  and the printing in OPAL-8 format, as well as some common attributes.
//
//  It defines a registry for attribute cells, used in the ATTLIST command
//  only.  The exemplar constructors for all OPAL element commands store all
//  defined attribute names in this registry.  The ATTLIST command can walk
//  through a beam line or sequence, and call the fillRegisteredAttributes()
//  method for each element.  This method will fill in the values for all
//  attributes which exist for this element, and the ATTLIST command can
//  look them up with findRegisteredAttribute() to build up a print line.

class OpalElement: public Element {

public:

  /// The common attributes for all elements.
  enum {
    TYPE,     // The design type.
    APERT,    // The aperture data.
    LENGTH,   // The element length.
    WAKEF,    // The wake function to be used
    COMMON
  };

  /// Switch for value desired on ATTLIST command.
  enum ValueFlag {
    ACTUAL_FLAG,      // Actual field values (design + error).
    IDEAL_FLAG,       // Ideal field values (design only).
    ERROR_FLAG        // Field errors.
  };

  virtual ~OpalElement();

  /// Read element from the DOOM data base.
  //  Access the element attributes via the DOOM index table.
  virtual void doomGet(const DoomReader &);

  /// Write element to the DOOM data base.
  //  Access the element attributes via the DOOM index table.
  virtual void doomPut(DoomWriter &) const;


  /// Fill in all registered attributes.
  virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

  /// Find a registered attribute.
  //  Return a pointer to the AttCell for a named attribute.
  static AttCell *findRegisteredAttribute(const string &name);


  /// Return element length.
  virtual double getLength() const;

  /// Return the element's type name.
  const string getTypeName() const;

  //return the element aperture vector
  std::vector<double> getApert() const;


  /// Return the element's type name.
  const string getWakeF() const;

  const string getWMaterial() const;

  const string getWakeGeom() const;

  std::vector<double> getWakeParam() const;

  const string getWakeConductivity() const;

  /// Parse the element.
  //  This special version for elements handles unknown attributes by
  //  appending them to the attribute list.
  virtual void parse(Statement &);

  /// Print the object.
  //  This special version handles special printing in OPAL-8 format.
  virtual void print(std::ostream &) const;

  /// Store a registered real attribute.
  static void setRegisteredAttribute(const string &, double);

  /// Store a registered string attribute.
  static void setRegisteredAttribute(const string &, const string &);

  /// Transmit the ``unknown'' (not known to OPAL) attributes to CLASSIC.
  virtual void updateUnknown(ElementBase *);

protected:

  /// Exemplar constructor.
  OpalElement(int size, const char *name, const char *help);

  /// Clone constructor.
  OpalElement(const string &name, OpalElement *parent);

  /// Print multipole components in OPAL-8 format.
  //  This function is accessible to all multipole-like elements
  //  (RBend, SBend, Quadrupole, Sextupole, Octupole, Multipole).
  static void printMultipoleStrength(std::ostream &os,
                                     int order,
                                     int &len,
                                     const string &sName,
                                     const string &tName,
                                     const Attribute &length,
                                     const Attribute &vNorm,
                                     const Attribute &vSkew);

  /// Print an attribute with a OPAL-8 name (as an expression).
  static void printAttribute(std::ostream &os,
                             const string &name,
                             const string &image,
                             int &len);

  /// Print an attribute with a OPAL-8 name (as a constant).
  static void printAttribute(std::ostream &os,
                             const string &name,
                             double value,
                             int &len);

  /// Register a ``real'' element attribute.
  //  A registered attribute can be listed by the ATTLIST command.
  static AttCell *registerRealAttribute(const string &name);

  /// Register a ``string'' element attribute.
  //  A registered attribute can be listed by the ATTLIST command.
  static AttCell *registerStringAttribute(const string &name);

  /// The registry for named attributes.
  static std::map < string, OwnPtr<AttCell> > attributeRegistry;

private:

  // Not implemented.
  OpalElement();
  void operator=(const OpalElement &);

  // The original size of the attribute list.
  int itsSize;
};

#endif // OPAL_OpalElement_HH
