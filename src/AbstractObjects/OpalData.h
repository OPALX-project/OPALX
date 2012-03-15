#ifndef OPAL_OpalData_HH
#define OPAL_OpalData_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalData.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalData
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------
#include "Structure/DataSink.h"
#include "AbstractObjects/ObjectFunction.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include <iosfwd>
#include <string>
#include <vector>
using std::string;

class AttributeBase;
class Object;
class Table;
class ValueDefinition;
class DataSink;

// Class OpalData
// ------------------------------------------------------------------------
/// The global OPAL structure.
//  The OPAL object holds all global data required for an OPAL execution.
//  In particular it contains the main Directory, which allows retrieval
//  of command objects by their name.  For other data refer to the
//  implementation file.

class OpalData {

public:

    /// Constructor.
    //  Initialise OPAL execution.
    OpalData();

    ~OpalData();

    /// Apply a function to all objects.
    //  Loop over the directory and apply the given functor object to each
    //  object in turn.
    void apply(const ObjectFunction &);

    /// Create new object.
    //  No replacement is allowed; if an object with the same name exists,
    //  throw [b]OpalException[/b].
    void create(Object *newObject);

    /// Define a new object.
    //  Replacement is allowed; however [b]OpalException[/b] is thrown,
    //  if the replacement cannot be done.
    void define(Object *newObject);

    /// Delete existing entry.
    //  Identified by [b]name[/b].
    void erase(const string &name);

    /// Find entry.
    //  Identified by [b]name[/b].
    Object *find(const string &name);

    /// Return value of global reference momentum.
    double getP0() const;

    /// Invalidate expressions.
    //  Force re-evaluation of all expressions before next command is
    //  executed.
    //  Also set the [b]modified[/b] flag in [b]object[/b], if not NULL.
    void makeDirty(Object *object);

    /// Print all objects.
    //  Loop over the directory and print each object whose name matches
    //  the regular expression [b]pattern[/b].
    void printNames(std::ostream &stream, const string &pattern);

    /// Register table.
    //  Register the table [b]t[/b].
    //  Registered tables are invalidated to be refilled when an object
    //  on which they depend is changed or replaced.
    void registerTable(Table *t);

    /// Unregister table.
    void unregisterTable(Table *t);

    /// Register expression.
    //  Registered expressions are invalidated to be recomputed when
    //  any object in the directory is changed or replaced.
    void registerExpression(AttributeBase *);

    /// Unregister expression.
    void unregisterExpression(AttributeBase *);

    /// Set the global momentum.
    void setP0(ValueDefinition *p0);

    /// Store the page title.
    void storeTitle(const string &);

    /// Print the page title.
    void printTitle(std::ostream &);

    /// Get the title string
    string getTitle();

    /// Update all objects.
    //  Loop over the directory and notify all objects to update themselves.
    void update();

    /// Clear Reference.
    //  This functor is used to clear the reference count stored in an object.
    struct ClearReference: public ObjectFunction {
        virtual void operator()(Object *) const;
    };

    std::vector<string> getAllNames();

    /// true if we do a restart run
    bool inRestartRun();

    /// set OPAL in restart mode
    void setRestartRun(const bool &value = true);

    /// store the location where to restart
    void setRestartStep(int s);

    /// get the step where to restart
    int getRestartStep();

    /// get opals input filename
    string getInputFn();

    /// store opals input filename
    void storeInputFn(const string &fn);

    /// get opals restart h5 format filename
    string getRestartFileName();

    /// store opals restart h5 format filename
    void setRestartFileName(string s);

    /// true if we do a restart from specified h5 file
    bool hasRestartFile();

    /// set the dump frequency as found in restart file
    void setRestartDumpFreq(const int &N);

    /// get the dump frequency as found in restart file
    int getRestartDumpFreq() const;

    /// set the last step in a run for possible follow-up run
    void setLastStep(const int &step);

    /// get the last step from a possible previous run
    int getLastStep() const;

    /// true if we already allocated a ParticleBunch object
    bool hasBunchAllocated();

    void bunchIsAllocated();

    PartBunch *getPartBunch();

    void setPartBunch(PartBunch *p);

    /// true if we already allocated a DataSink object
    bool hasDataSinkAllocated();

    DataSink *getDataSink();

    void setDataSink(DataSink *s);

    /// true if we already allocated a ParticleBunch object
    bool hasSLBunchAllocated();

    void slbunchIsAllocated();

    EnvelopeBunch *getSLPartBunch();

    void setSLPartBunch(EnvelopeBunch *p);

private:

    // Not implemented.
    OpalData(const OpalData &);
    void operator=(const OpalData &);

    // The private implementation details.
    struct OpalDataImpl *p;
};

extern OpalData OPAL;

#endif // OPAL_OpalData_HH
