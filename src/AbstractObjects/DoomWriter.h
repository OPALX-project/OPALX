#ifndef OPAL_DoomWriter_HH
#define OPAL_DoomWriter_HH

// ------------------------------------------------------------------------
// $RCSfile: DoomWriter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomWriter
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


// Class DoomWriter
// ------------------------------------------------------------------------
/// Encapsulates writing to the DOOM data base.
//  This class holds a DOOM ``object'' which contains several names, an
//  array of real values, an array of integers, and an array of C strings.
//  It is constructed by calling member functions of DoomWriter.  The
//  destructor of DoomWriter writes this object automatically to the DOOM
//  data base when the DoomWriter goes out of scope.

class DoomWriter {

public:

    /// Default constructor.
    //  Constructs an empty writer object.
    //  An object constructed with the default constructor can only be
    //  written if its name is set before the destructor is called.
    DoomWriter();

    /// Constructor.
    //  Construct an empty writer object for the OPAL Object [b]name[/b].
    explicit DoomWriter(const string &name);

    /// Constructor.
    //  Construct an empty writer object for the OPAL Object [b]name[/b].
    //  Like [b]DoomWriter(name)[/b], but use [b]keyList[/b] to further
    //  identify the object.
    DoomWriter(const string &name, const int keyList[]);

    /// Destructor.
    //  Save the object to the DOOM data base, if its name is set.
    ~DoomWriter();

    /// Store real value.
    //  Identified by [b]index[/b].
    void putReal(int index, double);

    /// Store real value.
    //  Identified by [b]index[/b].
    void putInt(int index, int);

    /// Store real value.
    //  Identified by [b]index[/b].
    void putString(int index, const string &);

    /// Set object name.
    //  This is the name of the OPAL object itself.
    void setObjectName(const string &name);

    /// Set object name.
    //  Like [b]setObjectName(name)[/b], but use [b]keyList[/b] to further
    //  identify the object.
    void setObjectName(const string &name, const int keyList[]);

    /// Set base name.
    //  This is the name of the OPAL exemplar object from which the object
    //  is ultimately derived.
    void setBaseName(const string &name);

    /// Set parent name.
    //  This is the name of the OPAL class object from which the object
    //  is derived directly.
    void setParentName(const string &name);

    /// Set type name.
    //  This is the DOOM object type like "ELEMENT", "ACTION", etc.
    void setTypeName(const string &name);

private:

    // Not implemented.
    DoomWriter(const DoomWriter &);
    void operator=(const DoomWriter &);

    // Truncate name to 24 characters.
    void truncateName(const string &name, char nam[24]);

    // The encapsulated DOOM data.
    struct object *doomStruct;
};

#endif // OPAL_DoomWriter_HH
