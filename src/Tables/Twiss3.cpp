// ------------------------------------------------------------------------
// $RCSfile: Twiss3.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Twiss3
//   The class for OPAL TWISS3 commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/Twiss3.h"
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

using std::setw;


// Class Twiss3
// ------------------------------------------------------------------------

// The attributes of class Twiss3.
namespace {
    enum {
        TABLE,       // The name of the table to be listed.
        FNAME,       // The name of the file to be written.
        SIZE
    };
}


Twiss3::Twiss3():
    Action(SIZE, "TWISS3",
           "The \"TWISS3\" statement lists a named \"TWISS\" table in "
           "Mais-Ripken representation.") {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of table to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "TWISS3");
}


Twiss3::Twiss3(const string &name, Twiss3 *parent):
    Action(name, parent)
{}


Twiss3::~Twiss3()
{}


Twiss3 *Twiss3::clone(const string &name) {
    return new Twiss3(name, this);
}


void Twiss3::execute() {
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
                throw OpalException("Twiss3::execute()",
                                    "Unable to open output stream \"" +
                                    fileName + "\".");
            }
        }
    } else {
        throw OpalException("Twiss3::execute()",
                            "Twiss table \"" + tableName + "\" not found.");
    }
}


void Twiss3::format(std::ostream &os, const Twiss *table) {
    if(Options::tfsFormat) {
        formatTFS(os, table);
    } else {
        formatPrint(os, table);
    }
}


void Twiss3::formatPrint(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Write table specific header.
    table->printTableTitle(os, "Mais-Ripken lattice functions");
    os << string(124, '-') << '\n';
    os << "Element" << string(24, ' ') << "S"
       << string(10, ' ') << "XC" << string(9, ' ') << "PXC"
       << string(10, ' ') << "YC" << string(9, ' ') << "PYC"
       << string(10, ' ') << "TC" << string(9, ' ') << "PTC\n";
    os << "Mode          MU"
       << "        BETX        GAMX        ALFX"
       << "        BETY        GAMY        ALFY"
       << "        BETT        GAMY        ALFT\n";
    os << string(124, '-') << '\n';
    // Jumbled function names affecting PRINT style output
    // fixed at 15:26:46 on 9 Aug 2000 by JMJ

    // Write table body.
    for(Twiss::TLine::const_iterator row = table->begin();
        row != table->end(); ++row) {
        if(row->getSelectionFlag()) {
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
            os << setw(16) << table->getS(*row);

            FVector<double, 6> orbit = table->getOrbit(*row);
            for(int i = 0; i < 6; ++i) {
                os << setw(12) << orbit[i];
            }
            os << '\n';

            for(int mode = 0; mode < 3; ++mode) {
                os << setw(4) << (mode + 1) << setw(12)
                   << table->getMUi(*row, mode);
                for(int plane = 0; plane < 3; ++plane) {
                    os << setw(12) << table->getBETik(*row, plane, mode)
                       << setw(12) << table->getGAMik(*row, plane, mode)
                       << setw(12) << table->getALFik(*row, plane, mode);
                }
                os << '\n';
            }
        }
    }

    os << string(124, '-') << std::endl;

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}


void Twiss3::formatTFS(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(12);
    std::ios::fmtflags old_flag = os.setf(std::ios::fixed, std::ios::floatfield);

    // Write table descriptors.
    os << "@ TYPE     %s  TWISS3\n";
    table->tfsTwissDescriptors(os);

    // Write column header names.
    std::string coNames[] = { " XC", " PXC", " YC", " PYC", " TC", " PTC" };
    os << "* NAME S";
    for(int i = 0; i < 6; ++i) os << coNames[i];
    static const char c[] = "XYT";
    for(int i = 0; i < 3; ++i) {
        for(int j = 1; j <= 3; ++j) {
            os << " BET" << c[i] << j << " ALF" << c[i] << j << " GAM" << c[i] << j;
        }
    }
    os << '\n';

    // Write column header formats.
    os << "$ %s";
    for(int i = 1; i <= 34; ++i) os << " %le";
    os << '\n';

    // Write table body.
    for(Twiss::TLine::const_iterator row = table->begin();
        row != table->end(); ++row) {
        if(row->getSelectionFlag()) {
            os << "  " << row->getElement()->getName()
               << ' ' << table->getS(*row);
            FVector<double, 6> orbit = table->getOrbit(*row);
            for(int i = 0; i < 6; ++i) {
                os << ' ' << orbit[i];
            }
            for(int i = 0; i < 3; ++i) {
                for(int j = 0; j < 3; ++j) {
                    os << ' ' << table->getBETik(*row, i, j)
                       << ' ' << table->getALFik(*row, i, j)
                       << ' ' << table->getGAMik(*row, i, j);
                }
            }
            os << '\n';
        }
    }

    // Restore the formatting flags.
    os.flush();
    os.precision(old_prec);
    os.flags(old_flag);
}
