// ------------------------------------------------------------------------
// $RCSfile: AttList.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttList
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:21 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Tables/AttList.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Table.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Elements/AttCell.h"
#include "Elements/OpalElement.h"
#include "Tables/AttWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include <fstream>
#include <iostream>
#include <vector>

using std::vector;


// Class AttList
// ------------------------------------------------------------------------

namespace {

    // The attributes of class AttList.
    enum {
        LINE,        // The name of the line to be listed.
        FNAME,       // The name of the file to be written.
        ALL,         // If true, list all columns.
        VALUE,       // Which value is desired: "ACTUAL", "IDEAL", "ERROR".
        COLUMN,      // The columns to be written.
        SIZE
    };
}


AttList::AttList():
    Action(SIZE, "ATTLIST",
           "The \"ATTLIST\" statement lists the element strengths "
           "in a beam line or a sequence.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of line to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "ATTLIST");
    itsAttr[ALL] = Attributes::makeBool
                   ("ALL", "Are all columns desired?");
    itsAttr[VALUE] = Attributes::makeString
                     ("VALUE", "Which value is desired: ACTUAL, IDEAL, or ERROR.", "ACTUAL");
    itsAttr[COLUMN] = Attributes::makeStringArray
                      ("COLUMN", "The columns to be written");
}


AttList::AttList(const string &name, AttList *parent):
    Action(name, parent)
{}


AttList::~AttList()
{}


AttList *AttList::clone(const string &name) {
    return new AttList(name, this);
}


void AttList::execute() {
    // Find beam sequence  or table definition.
    const string name = Attributes::getString(itsAttr[LINE]);
    const Beamline *line = 0;

    if(Object *obj = OpalData::getInstance()->find(name)) {
        if(BeamSequence *beamLine = dynamic_cast<BeamSequence *>(obj)) {
            line = beamLine->fetchLine();
        } else if(Table *table = dynamic_cast<Table *>(obj)) {
            line = table->getLine();
        } else {
            throw OpalException("Select::execute()",
                                "You cannot do an \"ATTLIST\" on \"" + name +
                                "\", it is neither a line nor a table.");
        }
    } else {
        throw OpalException("Select::execute()",
                            "Object \"" + name + "\" not found.");
    }

    // Select the file to be used.
    const string &fileName = Attributes::getString(itsAttr[FNAME]);
    if(fileName == "TERM") {
        writeTable(*line, std::cout);
    } else {
        std::ofstream os(fileName.c_str());

        if(os.good()) {
            writeTable(*line, os);
        } else {
            throw OpalException("AttList::execute()",
                                "Unable to open output stream \"" +
                                fileName + "\".");
        }
    }
}


void AttList::writeTable(const Beamline &line, std::ostream &os) {
    // Type of values desired.
    const string &value = Attributes::getString(itsAttr[VALUE]);
    OpalElement::ValueFlag flag = OpalElement::ACTUAL_FLAG;
    if(value == "ACTUAL") {
        flag = OpalElement::ACTUAL_FLAG;
    } else if(value == "IDEAL") {
        flag = OpalElement::IDEAL_FLAG;
    } else if(value == "ERROR") {
        flag = OpalElement::ERROR_FLAG;
    } else {
        throw OpalException("AttList::writeTable()",
                            "Unknown \"VALUE\" type \"" + value + "\".");
    }

    // Construct column access table.
    // This may throw, if a column is unknown.
    vector<string> header = Attributes::getStringArray(itsAttr[COLUMN]);
    vector<string>::size_type n = header.size();
    vector<AttCell *> buffer(n);
    for(vector<string>::size_type i = 0; i < n; ++i) {
        buffer[i] = OpalElement::findRegisteredAttribute(header[i]);
    }

    // Write table descriptors.
    OPALTimer::Timer timer;
    os << "@ TYPE     %s  ATTRIBUTE\n"
       << "@ LINE     %s  " << line.getName() << "\n"
       << "@ DATE     %s  " << timer.date() << "\n"
       << "@ TIME     %s  " << timer.time() << "\n"
       << "@ ORIGIN   %s  OPAL_9.5/4\n"
       << "@ COMMENT  %s  \""
       << "@ VALUE    %s  " << value << "\n";
    OpalData::getInstance()->printTitle(os);
    os << "\"\n";

    // Write column header names.
    os << '*';
    for(vector<string>::size_type i = 0; i < n; ++i) {
        os << ' ' << header[i];
    }
    os << '\n';

    // Write column header formats.
    os << '$';
    for(vector<string>::size_type i = 0; i < n; ++i) {
        os << ' ';
        buffer[i]->printFormat(os);
        buffer[i]->clearValue();
    }
    os << '\n';

    // List the table body.
    AttWriter writer(line, os, flag, buffer);
    writer.execute();
}
