// ------------------------------------------------------------------------
// $RCSfile: DoomDB.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomDB
//   Encapsulates the DOOM data base.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Fields/BMultipoleField.h"
#include "Utilities/OpalException.h"
#include "doomex.h"
#include <iostream>

using std::cerr;
using std::endl;


// Class DoomDB::Writer
// ------------------------------------------------------------------------

void DoomDB::Writer::operator()(Object *object) const {
    if(! object->isBuiltin()) {
        itsWriter.putString(size++, object->getOpalName());
        if(object->isDirty()) {
            DOOM_DB.writeObject(object);
        }
    }
}


// Class DoomDB
// ------------------------------------------------------------------------

DoomDB::DoomDB():
    dataBaseName(), dataBaseOpen(false),
    update(true), debugFlag(0), doomTime(0.0)
{}


DoomDB::~DoomDB()
{}


int DoomDB::getAttributeIndex(const string &name) {
    int low  = 0;
    int high = el_att_num - 1;

    while(low <= high) {
        int mid = (low + high) / 2;
        int pos = strcmp(name.c_str(), el_att_names[mid]);

        if(pos < 0) {
            high = mid - 1;
        } else if(pos > 0) {
            low  = mid + 1;
        } else {
            return el_att_pos[mid];
        }
    }

    return -1;
}


void DoomDB::firstOpen(const char *dbName) {
    dataBaseName = dbName;
    open();
}


void DoomDB::reOpen() {
    open();
}


void DoomDB::shut() {
    close();
}


void DoomDB::setUpdate(bool flag) {
    update = flag;
}


void DoomDB::setDebug(int flag) {
    debugFlag = flag;
}


Object *DoomDB::readObject(const string &name) {
    Object *obj = OPAL.find(name);

    if(dataBaseOpen) {
        char key[24];
        truncateName(key, name);
        double DBtime;
        int zero = 0;
        doom_gtime(key, &zero, &DBtime);

        if(obj == 0  ||  DBtime > obj->getDoomTime()) {
            if(debugFlag > 0) {
                cerr << endl << "reading: " << name << endl << endl;
            }

            // Read the DOOM data for this object.
            DoomReader reader(name);

            // Make sure the parent object is available.
            Object *parent = readObject(reader.getParentName());

            // Create clone of parent and extract its data.
            obj = parent->clone(name);
            obj->doomGet(reader);
            obj->setDoomTime(DBtime);
            OPAL.define(obj);

            // Reset the dirty flag which was set by OPAL.define(obj).
            obj->setDirty(false);
        }
    }

    return obj;
}


void DoomDB::writeObject(Object *object) {
    if(dataBaseOpen && object->isDirty()) {
        // Make sure parent is an exemplar, or exists in data base.
        Object *parent = object->getParent();
        writeObject(parent);

        // Write Object.
        if(debugFlag > 0) {
            cerr << endl << "writing: " << object->getOpalName() << endl << endl;
        }

        {
            // Save object to data base.
            DoomWriter writer(object->getOpalName());
            writer.setParentName(parent->getOpalName());
            writer.setBaseName(object->getBaseObject()->getOpalName());
            writer.setTypeName(object->getCategory());
            object->doomPut(writer);
            // The actual save happens in the destructor; requires inner scope.
        }

        // Reset the dirty flag, since the object is now the same in the data base.
        object->setDirty(false);
    }
}


bool DoomDB::readAlign(const string &name, int occur, Euclid3D &offset) {
    // Set up object key.
    char key[24];
    truncateName(key, name);

    // Fetch the errors from DOOM data base.
    int errors = 6;
    double align[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    doom_galign(key, &occur, &errors, align);

    // Construct the misalignment.
    if(errors > 0) {
        offset = Euclid3D(align[0], align[1], align[2],
                          align[3], align[4], align[5]);
        return true;
    } else {
        return false;
    }
}


void DoomDB::writeAlign(const string &name, int occur, const Euclid3D &offset) {
    // Set up object key.
    char key[24];
    truncateName(key, name);

    // Set up the misalign object.
    double align[6];
    offset.getAll(align[0], align[1], align[2], align[3], align[4], align[5]);

    // Save the errors to DOOM data base.
    int six = 6;
    doom_palign(key, &occur, &six, align);
}


bool DoomDB::readField(const string &name, int occur, double length,
                       const BMultipoleField &,
                       BMultipoleField &errorField) {
    double errors[256];
    int size = sizeof(errors) / sizeof(double);
    char key[24];
    truncateName(key, name);
    doom_gfield(key, &occur, &size, errors);
    if(length == 0.0) length = 1.0;

    if(size > 0) {
        int top = (size + 1) / 2;

        for(int comp = 1; comp <= top; ++comp) {
            errorField.setNormalComponent(comp, errors[2*comp-2] / length);
            errorField.setSkewComponent(comp, errors[2*comp-1] / length);
        }

        return true;
    } else {
        return false;
    }
}


void DoomDB::writeField(const string &name, int occur, double length,
                        const BMultipoleField &,
                        const BMultipoleField &errorField) {
    double errors[256];
    int size = 2 * errorField.order();
    if(length == 0.0) length = 1.0;

    for(int comp = 1; comp <= errorField.order(); ++comp) {
        errors[2*comp-2] = length * errorField.getNormalComponent(comp);
        errors[2*comp-1] = length * errorField.getSkewComponent(comp);
    }

    char key[24];
    truncateName(key, name);
    doom_pfield(key, &occur, &size, errors);
}


void DoomDB::setEnvironment(const char name[], const string &value) {
    char t1[24], t2[24];
    strcpy(t1, name);
    truncateName(t2, value);
    doom_setenv(t1, t2);
}


void DoomDB::close() {
    if(dataBaseOpen) {
        // Make sure all parameters are updated.
        OPAL.update();

        // Save modified directory.
        {
            DoomWriter directory("DIRECTORY");
            directory.setParentName("DIRECTORY");
            directory.setTypeName("DIRECTORY");
            DoomDB::Writer writer(directory);
            OPAL.apply(writer);
            // Directory is saved in destructor.
        }

        // Now the directory is saved; close the data base.
        doom_close();
        doom_time(&doomTime);

        if(debugFlag > 0) {
            cerr << endl << "Data base closed" << endl << endl;
        }

        dataBaseOpen = false;
    }
}


void DoomDB::open() {
    if(! dataBaseOpen  &&  ! dataBaseName.empty()) {
        try {
            // Load complete directory.
            // The const_cast<> is required, but harmless.
            doom_open(const_cast<char *>(dataBaseName.c_str()));
            dataBaseOpen = true;
            int zero = 0;
            double DBtime;
            doom_gtime("DIRECTORY", &zero, &DBtime);

            if(DBtime > doomTime) {
                DoomReader directory("DIRECTORY");
                int top = directory.getStringSize();

                for(int i = 0; i < top; i++) {
                    readObject(directory.getString(i));
                }
            }

            if(debugFlag > 0) {
                cerr << endl << "Data base opened." << endl << endl;
            }

            dataBaseOpen = true;
        } catch(OpalException &ex) {
            cerr << endl << "*** Error detected by function \""
                 << ex.where() << "\"" << endl;
            cerr << "    " << ex.what() << endl << endl;
        }
    }
}


void DoomDB::truncateName(char target[], const string &source) {
    int length = source.length();
    if(length >= 24) length = 23;
    strncpy(target, source.data(), length);
    target[length] = '\0';
}
