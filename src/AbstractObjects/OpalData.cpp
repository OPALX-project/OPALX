//
// Class OpalData
//   The global OPAL structure.
//   The OPAL object holds all global data required for a OPAL execution.
//   In particular it contains the main Directory, which allows retrieval
//   of command objects by their name.  For other data refer to the
//   implementation file.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbstractObjects/OpalData.h"

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Directory.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/ObjectFunction.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Attributes/Attributes.h"
#include "OpalParser/FileStream.h"
#include "OpalParser/OpalParser.h"
#include "OpalParser/StringStream.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Structure/DataSink.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <algorithm>
#include <iostream>
#include <list>
#include <set>

// Struct OpalDataImpl.
// ------------------------------------------------------------------------
struct OpalDataImpl {
    OpalDataImpl();
    ~OpalDataImpl();

    // The main object directory.
    Directory mainDirectory;

    // The keywords, by name, with the function that makes an object of each.
    std::map<std::string, OpalData::Maker> keywords;

    // The value of the global momentum.
    ValueDefinition* referenceMomentum;

    // The flag telling that something has changed.
    bool modified;

    // The set of expressions to be invalidated when something changes.
    std::set<AttributeBase*> exprDirectory;
    typedef std::set<AttributeBase*>::iterator exprIterator;

    // true if we restart a simulation
    bool isRestart_m;

    // Where to resume in a restart run
    std::string restartFn_m;

    // Input file name
    std::string inputFn_m;

    std::set<std::string> outFiles_m;

    /// Mode for writing files
    OpalData::OpenMode openMode_m = OpalData::OpenMode::WRITE;

    bool hasDataSinkAllocated_m;

    DataSink* dataSink_m;

    std::vector<MaxPhasesT> maxPhases_m;

    bool isInPrepState_m;

    std::vector<std::string> arguments_m;
};

OpalDataImpl::OpalDataImpl()
    : mainDirectory(),
      referenceMomentum(0),
      modified(false),
      isRestart_m(false),
      hasDataSinkAllocated_m(false),
      isInPrepState_m(false) {
    dataSink_m = nullptr;
}

OpalDataImpl::~OpalDataImpl() {
    // Make sure the main directory is cleared before the expressions are deleted.
    delete dataSink_m;

    mainDirectory.erase();
    exprDirectory.clear();
}

bool OpalData::isInstantiated = false;
OpalData* OpalData::instance  = nullptr;

OpalData* OpalData::getInstance() {
    if (!isInstantiated) {
        instance       = new OpalData();
        isInstantiated = true;
        return instance;
    } else {
        return instance;
    }
}

void OpalData::deleteInstance() {
    delete instance;
    instance       = nullptr;
    isInstantiated = false;
}

OpalData::OpalData() { p = new OpalDataImpl(); }

OpalData::~OpalData() { delete p; }

bool OpalData::isInPrepState() { return p->isInPrepState_m; }

void OpalData::setInPrepState(bool state) { p->isInPrepState_m = state; }

OpalData::PreparationState::PreparationState(OpalData& opalData) : opalData_m(opalData) {
    opalData_m.setInPrepState(true);
}

OpalData::PreparationState::~PreparationState() noexcept { finish(); }

void OpalData::PreparationState::finish() noexcept {
    if (active_m) {
        opalData_m.setInPrepState(false);
        active_m = false;
    }
}

OpalData::PreparationState OpalData::enterPreparationState() { return PreparationState(*this); }

bool OpalData::inRestartRun() { return p->isRestart_m; }

void OpalData::setRestartRun(const bool& value) { p->isRestart_m = value; }

std::string OpalData::getRestartFileName() { return p->restartFn_m; }

void OpalData::setRestartFileName(std::string s) { p->restartFn_m = s; }

void OpalData::setOpenMode(OpenMode openMode) { p->openMode_m = openMode; }

OpalData::OpenMode OpalData::getOpenMode() const { return p->openMode_m; }

bool OpalData::hasDataSinkAllocated() { return p->hasDataSinkAllocated_m; }

void OpalData::setDataSink(DataSink* s) {
    p->dataSink_m             = s;
    p->hasDataSinkAllocated_m = true;
}

DataSink* OpalData::getDataSink() { return p->dataSink_m; }

void OpalData::setMaxPhase(std::string elName, double phi) {
    p->maxPhases_m.push_back(MaxPhasesT(elName, phi));
}

std::vector<MaxPhasesT>::iterator OpalData::getFirstMaxPhases() { return p->maxPhases_m.begin(); }

std::vector<MaxPhasesT>::iterator OpalData::getLastMaxPhases() { return p->maxPhases_m.end(); }

int OpalData::getNumberOfMaxPhases() { return p->maxPhases_m.size(); }

void OpalData::apply(const ObjectFunction& fun) {
    for (ObjectDir::iterator i = p->mainDirectory.begin(); i != p->mainDirectory.end(); ++i) {
        fun(&*i->second);
    }
}

void OpalData::create(Object* newObject) {
    // Test for existing node with same name.
    const std::string name = newObject->getOpalName();
    Object* oldObject      = p->mainDirectory.find(name);

    if (oldObject != nullptr || isKeyword(name)) {
        throw OpalException(
                "OpalData::create()", "You cannot replace the object \"" + name + "\".");
    } else {
        p->mainDirectory.insert(name, newObject);
    }
}

void OpalData::define(Object* newObject) {
    // Test for existing node with same name.
    const std::string name = newObject->getOpalName();
    Object* oldObject      = p->mainDirectory.find(name);

    // A keyword cannot be redefined.
    if (isKeyword(name)) {
        throw OpalException(
                "OpalData::define()", "You cannot replace the object \"" + name + "\".");
    }

    if (oldObject != nullptr && oldObject != newObject) {
        // Attempt to replace an object.
        if (oldObject->isBuiltin() || !oldObject->canReplaceBy(newObject)) {
            throw OpalException(
                    "OpalData::define()", "You cannot replace the object \"" + name + "\".");
        } else {
            if (Options::info) {
                *ippl::Info << "Replacing the object \"" << name << "\"." << endl;
            }

            // Replace all references to this object.
            for (ObjectDir::iterator i = p->mainDirectory.begin(); i != p->mainDirectory.end();
                 ++i) {
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
    if (name == "P0") {
        if (ValueDefinition* p0 = dynamic_cast<ValueDefinition*>(newObject)) {
            setP0(p0);
        }
    }
}

void OpalData::erase(const std::string& name) {
    Object* oldObject = p->mainDirectory.find(name);

    if (oldObject != nullptr) p->mainDirectory.erase(name);
}

Object* OpalData::find(const std::string& name) { return p->mainDirectory.find(name); }

void OpalData::addKeyword(const std::string& keyword, Maker maker) {
    if (p->keywords.count(keyword) != 0 || p->mainDirectory.find(keyword) != nullptr) {
        throw OpalException(
                "OpalData::addKeyword()", "You cannot replace the object \"" + keyword + "\".");
    }
    p->keywords.emplace(keyword, std::move(maker));
}

bool OpalData::isKeyword(const std::string& name) const { return p->keywords.count(name) != 0; }

Object* OpalData::makeObject(const std::string& keyword, const std::string& name) {
    const auto entry = p->keywords.find(keyword);
    if (entry == p->keywords.end()) return nullptr;
    Object* object = entry->second();
    object->setOpalName(name);
    return object;
}

Object* OpalData::copyObject(const Object& source, const std::string& name) {
    Object* copy = makeObject(source.getKeyword(), name);
    if (copy == nullptr) {
        throw OpalException(
                "OpalData::copyObject()", "\"" + source.getOpalName() + "\" cannot be copied.");
    }
    copy->copyAttributes(source);
    return copy;
}

double OpalData::getP0() const { return p->referenceMomentum->getReal() * Units::GeV2eV; }

void OpalData::makeDirty(Object* obj) {
    p->modified = true;
    if (obj) obj->setDirty(true);
}

void OpalData::registerExpression(AttributeBase* expr) { p->exprDirectory.insert(expr); }

void OpalData::unregisterExpression(AttributeBase* expr) { p->exprDirectory.erase(expr); }

void OpalData::setP0(ValueDefinition* p0) { p->referenceMomentum = p0; }

void OpalData::storeInputFn(const std::string& fn) { p->inputFn_m = fn; }

std::string OpalData::getAuxiliaryOutputDirectory() const { return "data"; }

std::string OpalData::getInputFn() { return p->inputFn_m; }

std::string OpalData::getInputBasename() {
    std::string& fn = p->inputFn_m;
    int const pdot  = fn.rfind(".");
    return fn.substr(0, pdot);
}

void OpalData::checkAndAddOutputFileName(const std::string& outfn) {
    if (p->outFiles_m.count(outfn) == 0) {
        p->outFiles_m.insert(outfn);
    } else {
        throw OpalException(
                "OpalData::checkAndAddOutputFileName",
                "Duplicate file name for output, '" + outfn + "', detected");
    }
}

void OpalData::update() {
    if (p->modified) {
        // Force re-evaluation of expressions.
        for (OpalDataImpl::exprIterator i = p->exprDirectory.begin(); i != p->exprDirectory.end();
             ++i) {
            (*i)->invalidate();
        }

        p->modified = false;
    }
}

void OpalData::storeArguments(int argc, char* argv[]) {
    p->arguments_m.clear();
    for (int i = 0; i < argc; ++i) {
        p->arguments_m.push_back(argv[i]);
    }
}

std::vector<std::string> OpalData::getArguments() { return p->arguments_m; }
