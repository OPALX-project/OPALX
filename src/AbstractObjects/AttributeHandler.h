#ifndef OPAL_AttributeHandler_HH
#define OPAL_AttributeHandler_HH

// ------------------------------------------------------------------------
// $RCSfile: AttributeHandler.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeHandler
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "MemoryManagement/RCObject.h"
#include "AbstractObjects/AttributeBase.h"
#include "MemoryManagement/Pointer.h"
#include <string>

using std::string;

class Attribute;
class AttributeBase;
class DoomReader;
class DoomWriter;
class Statement;


// Class AttributeHandler
// ------------------------------------------------------------------------
/// Abstract base class for attribute parsers.
//  An attribute parser defines the data type for an attribute.  It is
//  used to parse the attribute, to read the attribute from the DOOM data
//  base, and to write it to that data base.  It contains the name and
//  help text for the attribute.  Optionally it may also contain a default
//  value for the attribute.
//  [p]
//  When ``is_readonly'' is true, the attribute cannot be redefined by
//  the user.
//  [p]
//  When ``is_deferred'' is true, the attribute must be re-evaluated
//  whenever it is used.  This is the case for random error values.
//  When ``is_deferred'' is false, any expression for the attribute is
//  cached.  It is re-evaluated only when any other definition has changed.

class AttributeHandler: public RCObject {

public:

  /// Constructor.
  //  Assigns the attribute name [b]name[/b] and the help text [b]help[/b],
  //  as well as a possible default value [b]def[/b]for the attribute.
  AttributeHandler(const string &name, const string &help, AttributeBase *def);

  virtual ~AttributeHandler();

  /// Make clone.
  //  Attribute handlers are always shared, so this method should never
  //  be called.  It exists only to fulfill the requirements of the class
  //  [b]Pointer[/b].
  virtual AttributeHandler *clone() const;

  /// Read attribute from the DOOM data base.
  //  Uses the DoomReader [b]r[/b] for the object being read,
  //  and the position [b]i[/b] within this reader.
  //  Called by [b]Attribute::doomGet()[/b]
  virtual void doomGet(Attribute &a, const DoomReader &r, int i) const = 0;

  /// Write the attribute a to the DOOM data base.
  //  Uses the DoomWriter [b]w[/b] for the object being read,
  //  and the position [b]i[/b] within this writer.
  //  Called by [b]Attribute::doomPut()[/b]
  virtual void doomPut(const Attribute &a, DoomWriter &w, int i) const = 0;

  /// Return default value.
  //  Return the default value stored in this parser.
  virtual AttributeBase *getDefault() const;

  /// Return help string.
  virtual const string &getHelp() const;

  /// Return attribute name.
  virtual const string &getName() const;

  /// Return attribute type.
  //  Return a string describing the attribute type
  //  ("logical", "real", etc.).
  virtual const string &getType() const = 0;

  /// Parse new value.
  //  Parse value from the statement [b]s[/b] and assign it to the
  //  attribute [b]a[/b].
  virtual void parse(Attribute &a, Statement &s, bool eval) const = 0;

  /// Parse component value.
  //  Parse value from the statement [b]s[/b] and assign it to the
  //  attribute [b]a[/b], indexed by [b]i[/b].
  //  The default version assumes that the value is scalar,
  //  and it throws [b]OpalException[/b].
  virtual void parseComponent
  (Attribute &a, Statement &s, bool eval, int i) const;

  /// Return defer flag.
  //  True, if any expression evaluation is to be deferred.
  //  See [b]Expressions::ADeferred[/b] and [b]Expressions::SDeferred[/b]
  //  for details.
  bool isDeferred() const;

  /// Set or reset defer flag.
  //  If the flag is set, expressions are evaluated only when the value
  //  is fetched.
  void setDeferred(bool);

  /// Return read-only flag.
  //  If [b]parse[/b] is called with this flag set,
  //  then [b]OpalException[/b] is thrown.
  bool isReadOnly() const;

  /// Set or reset read-only flag.
  //  If [b]parse[/b] is called with the flag set,
  //  then [b]OpalException[/b] is thrown.
  void setReadOnly(bool);

protected:

  /// Attribute name.
  const string itsName;

  /// Help text.
  const string itsHelp;

  /// Default value.
  Pointer <AttributeBase> itsDefault;

  /// Defer flag.
  bool is_deferred;

  /// Read-only flag.
  bool is_readonly;

private:

  // Not implemented.
  AttributeHandler();
  AttributeHandler(const AttributeHandler &);
  void operator=(const AttributeHandler &);
};

#endif // OPAL_AttributeHandler_HH
