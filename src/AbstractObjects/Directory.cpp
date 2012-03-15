// ------------------------------------------------------------------------
// $RCSfile: Directory.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Directory
//   A directory for OPAL objects.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Directory.h"
#include "AbstractObjects/Object.h"


// Class Directory
// ------------------------------------------------------------------------

Directory::Directory():
    dir()
{}


Directory::~Directory() {
    erase();
}


ObjectDir::iterator Directory::begin() {
    return dir.begin();
}


ObjectDir::const_iterator Directory::begin() const {
    return dir.begin();
}


ObjectDir::iterator Directory::end() {
    return dir.end();
}


ObjectDir::const_iterator Directory::end() const {
    return dir.end();
}


void Directory::erase() {
    dir.erase(dir.begin(), dir.end());
}


void Directory::erase(const string &name) {
    dir.erase(name);
}


Object *Directory::find(const string &name) const {
    ObjectDir::const_iterator index = dir.find(name);

    if(index == dir.end()) {
        return 0;
    } else {
        return &*index->second;
    }
}


void Directory::insert(const string &name, Object *newObject) {
    ObjectDir::value_type p(name, Pointer<Object>(newObject));
    dir.insert(p);
}
