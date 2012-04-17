// ------------------------------------------------------------------------
// $RCSfile: OpalData.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalData
//   The global OPAL structure.
//   The OPAL object holds all global data required for a OPAL execution.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/OpalData.h"
#include "Structure/DataSink.h"
#include "Structure/BoundaryGeometry.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Directory.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/ObjectFunction.h"
#include "AbstractObjects/Table.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/RegularExpression.h"
#include <iostream>
#include <list>
#include <set>

// DTA
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/StringConstant.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Parser/StringStream.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/bet/EnvelopeBunch.h"
// /DTA

using namespace std;

// Class OpalData::ClearReference
// ------------------------------------------------------------------------

void OpalData::ClearReference::operator()(Object *object) const {
    object->clear();
}


// Struct OpalDataImpl.
// ------------------------------------------------------------------------

struct OpalDataImpl {
    OpalDataImpl();
    ~OpalDataImpl();

    // The main object directory.
    Directory mainDirectory;

    // The value of the global momentum.
    ValueDefinition *referenceMomentum;

    // The flag telling that something has changed.
    bool modified;

    // Directory of tables, for recalculation when something changes.
    std::list<Table *> tableDirectory;
    typedef std::list<Table *>::iterator tableIterator;

    // The set of expressions to be invalidated when something changes.
    std::set <AttributeBase *> exprDirectory;
    typedef std::set<AttributeBase *>::iterator exprIterator;

    // The page title from the latest TITLE command.
    string itsTitle;

    // true if we restart a simulation
    bool isRestart;

    // Input file name
    string inputFn;

    // Where to resume in a restart run
    int restartStep;

    // Where to resume in a restart run
    string restartFn;

    // true if the name of a restartFile is specified
    bool hasRestartFile_m;

    // dump frequency as found in restart file
    int restart_dump_freq_m;

    // last step of a run
    int last_step_m;

    bool hasBunchAllocated_m;

    bool hasDataSinkAllocated_m;


    // The particle bunch to be tracked.
    PartBunch *bunch_m;

    DataSink *dataSink_m;

    bool hasSLBunchAllocated_m;

    // The particle bunch to be tracked.
    EnvelopeBunch *slbunch_m;

    // In units of seconds. This is half of tEmission
    double gPhaseShift_m;

    BoundaryGeometry *bg_m;

    vector<MaxPhasesT> maxPhases_m;

    // The cartesian mesh
    Mesh_t *mesh_m;

    // The field layout f
    FieldLayout_t *FL_m;

    // The particle layout
    Layout_t *PL_m;
};


OpalDataImpl::OpalDataImpl():
    mainDirectory(), referenceMomentum(0), modified(false), itsTitle(),
    restart_dump_freq_m(1), last_step_m(0) {
    bunch_m = 0;
    slbunch_m = 0;
    dataSink_m = 0;
    bg_m = 0;
    mesh_m = 0;
    FL_m = 0;
    PL_m = 0;
}


OpalDataImpl::~OpalDataImpl() {
    // Make sure the main directory is cleared before the directories
    // for tables and expressions are deleted.
    delete mesh_m;// somehow this needs to be deleted first
	delete FL_m;
    //delete PL_m; //this gets deleted by FL_m

    delete bunch_m;
    delete slbunch_m;
    delete bg_m;
    delete dataSink_m;


    mainDirectory.erase();
    tableDirectory.clear();
    exprDirectory.clear();
}


// Class OpalData
// ------------------------------------------------------------------------

bool OpalData::isInstatiated = false;
OpalData *OpalData::instance = 0;

OpalData *OpalData::getInstance() {
    if(!isInstatiated) {
        instance = new OpalData();
        isInstatiated = true;
        return instance;
    } else {
        return instance;
    }
}

void OpalData::deleteInstance() {
    delete instance;
    instance = 0;
    isInstatiated = false;
}


OpalData::OpalData() {
    p = new OpalDataImpl();
    p->isRestart = false;
    p->hasRestartFile_m = false;
    p->hasBunchAllocated_m = false;
    p->hasDataSinkAllocated_m = false;
    p->hasSLBunchAllocated_m = false;
    p->gPhaseShift_m = 0.0;
    p->maxPhases_m.clear();
}

OpalData::~OpalData() {
    delete p;
}

void OpalData::reset() {
    delete p;
    p = new OpalDataImpl();
    p->isRestart = false;
    p->hasRestartFile_m = false;
    p->hasBunchAllocated_m = false;
    p->hasDataSinkAllocated_m = false;
    p->hasSLBunchAllocated_m = false;
    p->gPhaseShift_m = 0.0;
    p->maxPhases_m.clear();
}

bool OpalData::inRestartRun() {
    return p->isRestart;
}

void OpalData::setRestartRun(const bool &value) {
    p->isRestart = value;
}


void OpalData::setRestartStep(int s) {
    p->restartStep = s;
}

int OpalData::getRestartStep() {
    return p->restartStep;
}


string OpalData::getRestartFileName() {
    return p->restartFn;
}


void OpalData::setRestartFileName(string s) {
    p->restartFn = s;
    p->hasRestartFile_m = true;
}

bool OpalData::hasRestartFile() {
    return p->hasRestartFile_m;

}

void OpalData::setRestartDumpFreq(const int &N) {
    p->restart_dump_freq_m = N;
}

int OpalData::getRestartDumpFreq() const {
    return p->restart_dump_freq_m;
}

void OpalData::setLastStep(const int &step) {
    p->last_step_m = step;
}

int OpalData::getLastStep() const {
    return p->last_step_m;
}

bool OpalData::hasSLBunchAllocated() {
    return p->hasSLBunchAllocated_m;
}

void OpalData::slbunchIsAllocated() {
    p->hasSLBunchAllocated_m = true;
}

void OpalData::setSLPartBunch(EnvelopeBunch *b) {
    p->slbunch_m = b;
}

EnvelopeBunch *OpalData::getSLPartBunch() {
    return p->slbunch_m;
}

bool OpalData::hasBunchAllocated() {
    return p->hasBunchAllocated_m;
}

void OpalData::bunchIsAllocated() {
    p->hasBunchAllocated_m = true;
}

void OpalData::setPartBunch(PartBunch *b) {
    p->bunch_m = b;
}

PartBunch *OpalData::getPartBunch() {
    return p->bunch_m;
}


bool OpalData::hasDataSinkAllocated() {
    return p->hasDataSinkAllocated_m;
}

void OpalData::setDataSink(DataSink *s) {
    p->dataSink_m = s;
    p->hasDataSinkAllocated_m = true;
}

DataSink *OpalData::getDataSink() {
    return p->dataSink_m;
}

void OpalData::setMaxPhase(string elName, double phi) {
    p->maxPhases_m.push_back(MaxPhasesT(elName, phi));
}

vector<MaxPhasesT>::iterator OpalData::getFirstMaxPhases() {
    return p->maxPhases_m.begin();
}

vector<MaxPhasesT>::iterator OpalData::getLastMaxPhases() {
    return p->maxPhases_m.end();
}

int OpalData::getNumberOfMaxPhases() {
    return p->maxPhases_m.size();
}


Mesh_t* OpalData::getMesh() {
	return p->mesh_m;
}

FieldLayout_t* OpalData::getFieldLayout() {
	return p->FL_m;
}

Layout_t* OpalData::getLayout() {
	return p->PL_m;
}

void OpalData::setMesh(Mesh_t *mesh) {
	p->mesh_m = mesh;
}

void OpalData::setFieldLayout(FieldLayout_t *fieldlayout) {
	p->FL_m = fieldlayout;
}

void OpalData::setLayout(Layout_t *layout) {
	p->PL_m = layout;
}

void OpalData::setGlobalPhaseShift(double shift) {
    /// units: (sec)
    p->gPhaseShift_m = shift;
}

double OpalData::getGlobalPhaseShift() {
    /// units: (sec)
    return p->gPhaseShift_m;
}

void   OpalData::setGlobalGeometry(BoundaryGeometry *bg) {
    p->bg_m = bg;
}

BoundaryGeometry *OpalData::getGlobalGeometry() {
    return p->bg_m;
}

void OpalData::apply(const ObjectFunction &fun) {
    for(ObjectDir::iterator i = p->mainDirectory.begin();
        i != p->mainDirectory.end(); ++i) {
        fun(&*i->second);
    }
}


void OpalData::create(Object *newObject) {
    // Test for existing node with same name.
    const string name = newObject->getOpalName();
    Object *oldObject = p->mainDirectory.find(name);

    if(oldObject != 0) {
        throw OpalException("OpalData::create()",
                            "You cannot replace the object \"" + name + "\".");
    } else {
        p->mainDirectory.insert(name, newObject);
    }
}


void OpalData::define(Object *newObject) {
    // Test for existing node with same name.
    const string name = newObject->getOpalName();
    Object *oldObject = p->mainDirectory.find(name);

    if(oldObject != 0  &&  oldObject != newObject) {
        // Attempt to replace an object.
        if(oldObject->isBuiltin()  ||  ! oldObject->canReplaceBy(newObject)) {
            throw OpalException("OpalData::define()",
                                "You cannot replace the object \"" + name + "\".");
        } else {
            if(Options::info) {
                INFOMSG("Replacing the object \"" << name << "\"." << endl);
            }

            // Erase all tables which depend on the new object.
            OpalDataImpl::tableIterator i = p->tableDirectory.begin();
            while(i != p->tableDirectory.end()) {
                // We must increment i before calling erase(name),
                // since erase(name) removes "this" from "tables".
                Table *table = *i++;
                const string &tableName = table->getOpalName();

                if(table->isDependent(name)) {
                    if(Options::info) {
                        cerr << endl << "Erasing dependent table \"" << tableName
                             << "\"." << endl << endl;
                    }

                    // Remove table from directory.
                    // This erases the table from the main directory,
                    // and its destructor unregisters it from the table directory.
                    erase(tableName);
                }
            }

            // Replace all references to this object.
            for(ObjectDir::iterator i = p->mainDirectory.begin();
                i != p->mainDirectory.end(); ++i) {
                (*i).second->replace(oldObject, newObject);
            }

            // Remove old object.
            erase(name);
        }
    }

    // Force re-evaluation of expressions.
    p->modified = true;
    newObject->setDirty(true);
    p->mainDirectory.insert(name, newObject);

    // If this is a new definition of "P0", insert its definition.
    if(name == "P0") {
        if(ValueDefinition *p0 = dynamic_cast<ValueDefinition *>(newObject)) {
            setP0(p0);
        }
    }
}


void OpalData::erase(const string &name) {
    Object *oldObject = p->mainDirectory.find(name);

    if(oldObject != 0) {
        // Relink all children of "this" to "this->getParent()".
        for(ObjectDir::iterator i = p->mainDirectory.begin();
            i != p->mainDirectory.end(); ++i) {
            Object *child = &*i->second;
            if(child->getParent() == oldObject) {
                child->setParent(oldObject->getParent());
            }
        }

        // Remove the object.
        p->mainDirectory.erase(name);
    }
}


Object *OpalData::find(const string &name) {
    return p->mainDirectory.find(name);
}


double OpalData::getP0() const {
    static const double energy_scale = 1.0e+9;
    return p->referenceMomentum->getReal() * energy_scale;
}


void OpalData::makeDirty(Object *obj) {
    p->modified = true;
    if(obj) obj->setDirty(true);
}


void OpalData::printNames(std::ostream &os, const string &pattern) {
    int column = 0;
    RegularExpression regex(pattern);
    os << endl << "Object names matching the pattern \""
       << pattern << "\":" << endl;

    for(ObjectDir::const_iterator index = p->mainDirectory.begin();
        index != p->mainDirectory.end(); index++) {
        const string name = (*index).first;

        if(! name.empty()  &&  regex.match(name)) {
            os << name;

            if(column < 80) {
                column += name.length();

                do {
                    os << ' ';
                    column++;
                } while((column % 20) != 0);
            } else {
                os << endl;
                column = 0;
            }
        }
    }

    if(column) os << endl;
    os << endl;
}


void OpalData::registerTable(Table *table) {
    p->tableDirectory.push_back(table);
}


void OpalData::unregisterTable(Table *table) {
    for(OpalDataImpl::tableIterator i = p->tableDirectory.begin();
        i != p->tableDirectory.end();) {
        OpalDataImpl::tableIterator j = i++;
        if(*j == table) p->tableDirectory.erase(j);
    }
}


void OpalData::registerExpression(AttributeBase *expr) {
    p->exprDirectory.insert(expr);
}


void OpalData::unregisterExpression(AttributeBase *expr) {
    p->exprDirectory.erase(expr);
}


void OpalData::setP0(ValueDefinition *p0) {
    p->referenceMomentum = p0;
}


void OpalData::storeTitle(const string &title) {
    p->itsTitle = title;
}

void OpalData::storeInputFn(const string &fn) {
    p->inputFn = fn;
}


void OpalData::printTitle(std::ostream &os) {
    os << p->itsTitle;
}

string OpalData::getTitle() {
    return p->itsTitle;
}

string OpalData::getInputFn() {
    return p->inputFn;
}


void OpalData::update() {
    Inform msg("OpalData ");
    if(p->modified) {
        // Force re-evaluation of expressions.
        for(OpalDataImpl::exprIterator i = p->exprDirectory.begin();
            i != p->exprDirectory.end(); ++i) {
            (*i)->invalidate();
        }

        // Force refilling of dynamic tables.
        for(OpalDataImpl::tableIterator i = p->tableDirectory.begin();
            i != p->tableDirectory.end(); ++i) {
            (*i)->invalidate();
        }

        // Update all definitions.
        for(ObjectDir::iterator i = p->mainDirectory.begin();
            i != p->mainDirectory.end(); ++i) {
            (*i).second->update();
        }

        // Definitions are up-to-date.
        p->modified = false;
    }
}

std::vector<string> OpalData::getAllNames() {
    std::vector<string> result;

    for(ObjectDir::const_iterator index = p->mainDirectory.begin();
        index != p->mainDirectory.end(); index++) {
        string tmpName = (*index).first;
        if(!tmpName.empty()) result.push_back(tmpName);
        //// DTA
        //if (!tmpName.empty()) {
        //  Object *tmpObject = OPAL.find(tmpName);
        //  std::cerr << tmpObject->getCategory() << "\t" << tmpName << "\t";
        //  string bi = (tmpObject->isBuiltin()) ? "BUILT-IN" : "";
        //  std::cerr << bi << std::endl;
        //}
        //// /DTA
    }

    // DTA
    std::cout << "\nUser-defined variables:\n";
    const OpalParser mp;
    // FileStream *is;
    // is = new FileStream("/home/research/dabell/projects/opal9/src/tmp/myExpr.txt");
    // StringStream *ss;
    // ss = new StringStream("myVar:=ALPHA+BETA;");
    for(ObjectDir::const_iterator index = p->mainDirectory.begin();
        index != p->mainDirectory.end(); index++) {
        string tmpName = (*index).first;
        if(!tmpName.empty()) {
            Object *tmpObject = OpalData::getInstance()->find(tmpName);
            if(!tmpObject || tmpObject->isBuiltin()) continue;
            if(tmpObject->getCategory() == "VARIABLE") {
                std::cout << tmpName;
                if(dynamic_cast<RealVariable *>(&*tmpObject)) {
                    RealVariable &variable = *dynamic_cast<RealVariable *>(OpalData::getInstance()->find(tmpName));
                    std::cout << "\te= " << variable.value().getBase().getImage();
                    std::cout << "\tr= " << variable.getReal();
                    std::cout << "\ta= " << variable.itsAttr[0];
                    //Attributes::setReal(variable.itsAttr[0],137.);
                    //std::cout << "\te= " << variable.value().getBase().getImage();
                    //std::cout << "\tx= " << variable.itsAttr[0];
                }
                std::cout << std::endl;
            }
            //    if(tmpName=="MYL"){
            //      std::cout << "myL noted" << std::endl;
            //      RealVariable& variable = *dynamic_cast<RealVariable*>(OPAL.find(tmpName));
            //      std::cout << tmpName;
            //      std::cout << "\te= " << variable.value().getBase().getImage();
            //      std::cout << "\tr= " << variable.getReal();
            //      std::cout << "\ta= " << variable.itsAttr[0] << std::endl;
            //      mp.run(is);
            //      std::cout << tmpName;
            //      std::cout << "\te= " << variable.value().getBase().getImage();
            //      std::cout << "\tr= " << variable.getReal();
            //      std::cout << "\ta= " << variable.itsAttr[0] << std::endl;
            //      mp.run(ss);
            //    }
        }
    }
    std::cout << std::endl;
    ////TxAttributeSet data(variableName);
    ////const RealVariable& var = *dynamic_cast<RealVariable*>(OPAL.find(variableName));
    ////data.appendString("expression",variable.value().getBase().getImage());
    ////data.appendParam("value",variable.getReal());
    //// /DTA

    return result;
}

