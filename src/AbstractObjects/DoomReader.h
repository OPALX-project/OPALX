#ifndef OPAL_DoomReader_HH
#define OPAL_DoomReader_HH

// ------------------------------------------------------------------------
// $RCSfile: DoomReader.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomReader
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <string>

using std::string;

struct object;


// Class DoomReader
// ------------------------------------------------------------------------
/// Encapsulates reading from the DOOM data base.
//  This class reads an object from the DOOM data base and gives access
//  to its data.
//
//  The constructor fetches a named ``object'' from the data base and
//  keeps a pointer to that object.  The object contains several names,
//  an array of real values, an array of integers, and an array of
//  C strings. These can all be retrieved by calling member functions
//  of DoomReader.

class DoomReader {

public:

  /// Constructor.
  //  Constructs a reader object for the OPAL Object [b]name[/b].  It
  //  reads the data from the DOOM data base and stores them internally.
  explicit DoomReader(const string &name);

  /// Constructor.
  //  Like [b]DoomReader(name)[/b], but use [b]keyList[/b] to further
  //  identify the object.
  DoomReader(const string &name, int keyList[]);

  ~DoomReader();

  /// Number of real values.
  //  Return number of real values stored in this reader.
  int getRealSize() const; 

  /// Number of integer values.
  //  Return number of integer values stored in this reader.
  int getIntSize() const;

  /// Number of string values.
  //  Return number of string values stored in this reader.
  int getStringSize() const;

  /// Return real value.
  //  Identified by [b]index[/b].
  double getReal(int index) const;

  /// Return integer value.
  //  Identified by [b]index[/b].
  int getInt(int index) const;

  /// Return string value.
  //  Identified by [b]index[/b].
  string getString(int index) const;

  /// Return base name.
  //  This is the name of the OPAL exemplar object from which the object
  //  is ultimately derived.
  string getBaseName() const;

  /// Return parent name.
  //  This is the name of the OPAL class object from which the object
  //  is derived directly.
  string getParentName() const;

  /// Return type name.
  //  This is the DOOM object type like "ELEMENT", "ACTION", etc.
  string getTypeName() const;

private:

  // Not implemented.
  DoomReader();
  DoomReader(const DoomReader &);
  void operator=(const DoomReader &);

  // Method to read object from data base.
  void fill(char *key);

  // The encapsulated DOOM data.
  struct object *doomStruct;
};

#endif // OPAL_DoomReader_HH
