#ifndef OPAL_DoomDB_HH
#define OPAL_DoomDB_HH

// ------------------------------------------------------------------------
// $RCSfile: DoomDB.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomDB
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/ObjectFunction.h"
#include <string>

using std::string;

class BMultipoleField;
class DoomWriter;
class Euclid3D;
class Object;


// Class DoomDB
// ------------------------------------------------------------------------
/// Encapsulates the DOOM data base.
//  Objects are read from the data base using a DoomReader object,
//  and written using a DoomWriter object.  These objects encapsulate
//  Hans Grote's ``object'' structure.  Their constructors and destructors
//  handle the necessary transmission from and to the data base.  For more
//  details refer to those classes.

class DoomDB {

public:

    /// Constructor.
    //  Initialise the data base, but do not open it.
    DoomDB();

    /// Destructor.
    //  Do nothing.
    ~DoomDB();

    /// Find DOOM index for element attribute.
    //  Find the internal index for the element attribute identified
    //  by [b]name[/b].
    static int getAttributeIndex(const string &name);

    /// Open data base (first time) identified by file name [b]dbName[/b].
    void firstOpen(const char *dbName);

    /// Re-open data base which was already open.
    void reOpen();

    /// Close the data base.
    void shut();

    /// Set the debug flag.
    //  See the DOOM documentation for details.
    void setDebug(int);

    /// Set update flag.
    //  Called when an object is created or modified, to force writing
    //  all modified objects to the data base when OPAL is shut down.
    void setUpdate(bool);

    /// Read object.
    //  Identified by the name [b]name[/b].
    Object *readObject(const string &name);

    /// Write object.
    void writeObject(Object *object);

    /// Read misalignment.
    bool readAlign(const string &name, int occur, Euclid3D &offset);

    /// Write misalignment.
    void writeAlign(const string &name, int occur, const Euclid3D &offset);

    // Read multipole errors.
    // The fields in DOOM are integral (K_n / n!) * ds,
    // in OPAL-9 they are just K_n / n!, unless the length is zero,
    // in which case they are the same as in DOOM.
    static bool readField(const string &name, int occur, double length,
                          const BMultipoleField &mainField,
                          BMultipoleField &errorField);

    // Write multipole errors.
    // The fields in DOOM are integral (K_n / n!) * ds,
    // in OPAL-9 they are just K_n / n!, unless the length is zero,
    // in which case they are the same as in DOOM.
    static void writeField(const string &name, int occur, double length,
                           const BMultipoleField &mainField,
                           const BMultipoleField &errorField);

    /// Set DOOM environment variable.
    //  Call doom_setenv(name, value).
    //  See the DOOM documentation for details.
    static void setEnvironment(const char name[], const string &value);

private:

    // Ancillary class Writer.
    // Functor object.  OPAL saves objects to the doom data base by looping
    // over the OPAL directory and calling DoomDB::writer for each changed
    // object.

    struct Writer: ObjectFunction {
        Writer(DoomWriter &writer): itsWriter(writer), size(0) {}
        virtual void operator()(Object *) const;
        DoomWriter &itsWriter;
        mutable int size;
    };

    // Not implemented.
    DoomDB(const DoomDB &);
    void operator=(const DoomDB &);

    // Open/close the data base.
    void close();
    void open();

    // Copy a C++ string to a fixed-length char array, with terminating zero.
    static void truncateName(char target[], const string &source);

    // The data base name.
    string dataBaseName;
    bool dataBaseOpen;

    // The update and debug flags.
    bool update;
    int debugFlag;

    // The modification time.
    double doomTime;
};


/// The DOOM data base object.
//  This object encapsulates the DOOM data base for a particular OPAL run.
extern DoomDB DOOM_DB;

#endif // OPAL_DoomDB_HH
