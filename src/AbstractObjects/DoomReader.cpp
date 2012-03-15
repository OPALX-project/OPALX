// ------------------------------------------------------------------------
// $RCSfile: DoomReader.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomReader
//   This class reads an object from the DOOM data base and gives access
//   to its data.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/DoomReader.h"
#include "Utilities/OpalException.h"
#include "doom.h"
#include "doomex.h"
#include <cstring>

// Class DoomReader
// ------------------------------------------------------------------------


DoomReader::DoomReader(const string &name) {
    char nam[24];
    int len = name.length();
    if(len >= 24) len = 23;
    strncpy(nam, name.data(), len);
    nam[len] = '\0';
    fill(nam);
}


DoomReader::DoomReader(const string &name, int keyList[]) {
    char key[48];
    make_key(name.c_str(), keyList, key);
    fill(key);
}


DoomReader::~DoomReader() {
    // Do not destroy doomStruct, it is owned by the DOOM data base.
}


int DoomReader::getInt(int index) const {
    return (index < doomStruct->c_int) ? doomStruct->a_int[index] : 0;
}


double DoomReader::getReal(int index) const {
    return (index < doomStruct->c_dble) ? doomStruct->a_dble[index] : 0.0;
}


string DoomReader::getString(int index) const {
    return (index < doomStruct->c_obj) ? doomStruct->names[index] : "";
}


int DoomReader::getIntSize() const {
    return doomStruct->c_int;
}


int DoomReader::getRealSize() const {
    return doomStruct->c_dble;
}


int DoomReader::getStringSize() const {
    return doomStruct->c_obj;
}


string DoomReader::getBaseName() const {
    return doomStruct->base_name;
}


string DoomReader::getParentName() const {
    return doomStruct->par_name;
}


string DoomReader::getTypeName() const {
    return doomStruct->obj_type;
}


void DoomReader::fill(char *key) {
    doomStruct = doom_fetch(key);

    if(doomStruct == 0) {
        throw OpalException("DoomReader::fill()", "Unable to find key \"" +
                            string(key) + "\" in the current data base.");
    }
}
