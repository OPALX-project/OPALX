// ------------------------------------------------------------------------
// $RCSfile: Eigen.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Eigen
//   The class for OPAL EIGEN commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/Eigen.h"
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


// Class Eigen
// ------------------------------------------------------------------------

// The attributes of class Eigen.
namespace {
    enum {
        TABLE,       // The name of the table to be listed.
        FNAME,       // The name of the file to be written.
        SIZE
    };
}


Eigen::Eigen():
    Action(SIZE, "EIGEN",
           "The \"EIGEN\" statement lists the eigenvectors for a named "
           "\"TWISS\" table.") {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of table to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "EIGEN");
}


Eigen::Eigen(const string &name, Eigen *parent):
    Action(name, parent)
{}


Eigen::~Eigen()
{}


Eigen *Eigen::clone(const string &name) {
    return new Eigen(name, this);
}


void Eigen::execute() {
    string tableName = Attributes::getString(itsAttr[TABLE]);
    Twiss *table = dynamic_cast<Twiss *>(OPAL.find(tableName));

    if(table) {
        string fileName = Attributes::getString(itsAttr[FNAME]);
        if(fileName == "TERM") {
            format(std::cout, table);
        } else {
            std::ofstream os(fileName.c_str());

            if(os.good()) {
                format(os, table);
            } else {
                throw OpalException("Eigen::execute()",
                                    "Unable to open output stream \"" +
                                    fileName + "\".");
            }
        }
    } else {
        throw OpalException("Eigen::execute()",
                            "Twiss table \"" + tableName + "\" not found.");
    }
}


void Eigen::format(std::ostream &os, const Twiss *table) {
    if(Options::tfsFormat) {
        formatTFS(os, table);
    } else {
        formatPrint(os, table);
    }
}


void Eigen::formatPrint(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Print table specific header.
    table->printTableTitle(os, "Eigenvectors");
    os << string(118, '-') << '\n';
    os << "Element" << string(24, ' ') << "S       Orbit |"
       << string(25, ' ') << "E i g e n v e c t o r s\n";
    os << string(118, '-') << '\n';

    // Print table body.
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

            FMatrix<double, 6, 6> eigen = table->getCurlyA(*row);
            FVector<double, 6> orbit = table->getOrbit(*row);
            for(int i = 0; i < 6; ++i) {
                if(i != 0) os << string(32, ' ');
                os << std::setw(12) << orbit[i] << " |";
                for(int j = 0; j < 6; ++j) {
                    os << std::setw(12)  << eigen[i][j];
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


void Eigen::formatTFS(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(12);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Write table descriptors.
    os << "@ TYPE     %s  EIGEN\n";
    table->tfsTwissDescriptors(os);

    // Write column header names.
    os << "* NAME S XC PXC YC PYC TC PTC MU1 MU2 MU3";
    for(int i = 1; i <= 6; ++i) {
        for(int j = 1; j <= 6; ++j) os << " E" << i << j;
    }
    os << '\n';

    // Write column header formats.
    os << "$ %s";
    for(int i = 1; i <= 46; ++i) os << " %le";
    os << '\n';

    // Write table body.
    for(Twiss::TLine::const_iterator row = table->begin();
        row != table->end(); ++row) {
        if(row->getSelectionFlag()) {
            os << "  " << row->getElement()->getName()
               << ' ' << table->getS(*row, 0, 0);
            const FVector<double, 6> orbit = table->getOrbit(*row);
            for(int i = 0; i < 6; ++i) {
                os << ' ' << orbit[i];
            }

            os << ' ' << table->getMUi(*row, 0, 0)
               << ' ' << table->getMUi(*row, 1, 0)
               << ' ' << table->getMUi(*row, 2, 0);

            const FMatrix<double, 6, 6> eigen = table->getCurlyA(*row);
            for(int i = 0; i < 6; ++i) {
                for(int j = 0; j < 6; ++j) {
                    os << ' ' << eigen[i][j];
                }
            }
            os << '\n';
        }
    }

    // Restore the formatting flags.
    os.flush();
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}
