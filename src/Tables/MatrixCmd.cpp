// ------------------------------------------------------------------------
// $RCSfile: MatrixCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatrixCmd
//   The class for OPAL MATRIX commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/MatrixCmd.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Tables/Twiss.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <fstream>
#include <iomanip>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif
#include <iostream>

// Class MatrixCmd
// ------------------------------------------------------------------------

// The attributes of class MatrixCmd.
namespace {
    enum {
        TABLE,       // The name of the table to be listed.
        FNAME,       // The name of the file to be written.
        SIZE
    };
}


MatrixCmd::MatrixCmd():
    Action(SIZE, "MATRIX",
           "The \"MATRIX\" statement lists the accumulated transfer matrix "
           " for a named \"TWISS\" table.") {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of table to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "MATRIX");
}


MatrixCmd::MatrixCmd(const string &name, MatrixCmd *parent):
    Action(name, parent)
{}


MatrixCmd::~MatrixCmd()
{}


MatrixCmd *MatrixCmd::clone(const string &name) {
    return new MatrixCmd(name, this);
}


void MatrixCmd::execute() {
    string tableName = Attributes::getString(itsAttr[TABLE]);
    Twiss *table = dynamic_cast<Twiss *>(OpalData::getInstance()->find(tableName));

    if(table) {
        string fileName = Attributes::getString(itsAttr[FNAME]);
        if(fileName == "TERM") {
            format(std::cout, table);
        } else {
            std::ofstream os(fileName.c_str());

            if(os.good()) {
                format(os, table);
            } else {
                throw OpalException("MatrixCmd::execute()",
                                    "Unable to open output stream \"" +
                                    fileName + "\".");
            }
        }
    } else {
        throw OpalException("MatrixCmd::execute()",
                            "Twiss table \"" + tableName + "\" not found.");
    }
}


void MatrixCmd::format(std::ostream &os, const Twiss *table) {
        formatPrint(os, table);
}


void MatrixCmd::formatPrint(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Write table specific header.
    table->printTableTitle(os, "Accumulated transfer matrix");
    os << string(118, '-') << '\n';
    os << "Element" << string(24, ' ') << "S        Kick |"
       << string(25, ' ') << "T r a n s f e r   M a t r i x\n";
    os << string(118, '-') << '\n';

    // Write table body.
    for(Twiss::TLine::const_iterator row = table->begin();
        row != table->end(); ++row) {
        if(row->getSelectionFlag()) {
            os << '\n';
            string name = row->getElement()->getName();
            if(int occur = row->getCounter()) {
#if defined(__GNUC__) && __GNUC__ < 3
                char buffer[128];
                std::ostrstream tos(buffer, 128);
#else
                std::ostringstream tos;
#endif
                tos << name << '[' << occur << ']' << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
                name = buffer;
#else
                name = tos.str();
#endif
            }

            if(name.length() > 16) {
                // Truncate the element name.
                os << string(name, 0, 13) << ".. ";
            } else {
                // Left adjust the element name.
                os << name << string(16 - name.length(), ' ');
            }
            os << std::setw(16) << table->getS(*row);

            FVector<double, 6> orbit = table->getOrbit(*row);
            FMatrix<double, 6, 6> matrix = table->getMatrix(*row);
            for(int i = 0; i < 6; ++i) {
                if(i != 0) os << string(32, ' ');
                os << std::setw(12) << orbit[i] << " |";
                for(int j = 0; j < 6; ++j) {
                    os << std::setw(12) << matrix[i][j];
                }
                os << '\n';
            }
        }
    }

    os << string(118, '-') << std::endl;

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}


