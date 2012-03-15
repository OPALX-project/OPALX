// ------------------------------------------------------------------------
// $RCSfile: Insertion.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Insertion
//   The class for the OPAL TWISSTRACK command.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/Insertion.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/TableRowRep.h"
#include "Algorithms/Mapper.h"
#include "Attributes/Attributes.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/LinearMap.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include <iomanip>
#include <iostream>

using std::setw;


// Class Insertion
// ------------------------------------------------------------------------

Insertion::Insertion():
    Twiss(SIZE, "TWISSTRACK",
          "The \"TWISSTRACK\" command defines a table of lattice functions\n"
          "which can be matched or tabulated for an insertion.") {
    itsAttr[INIT] = Attributes::makeTableRow
                    ("INIT", "If given, this table position is used to initialise\n"
                     "\t\t\t\tif not given, use the following parameters to initialise");
    itsAttr[BETX] = Attributes::makeReal
                    ("BETX", "Initial horizontal beta", 1.0);
    itsAttr[ALFX] = Attributes::makeReal
                    ("ALFX", "Initial horizontal alpha");
    itsAttr[BETY] = Attributes::makeReal
                    ("BETY", "Initial vertical beta", 1.0);
    itsAttr[ALFY] = Attributes::makeReal
                    ("ALFY", "Initial vertical alpha");
    itsAttr[DX] = Attributes::makeReal
                  ("DX", "Initial horizontal dispersion");
    itsAttr[DPX] = Attributes::makeReal
                   ("DPX", "Initial slope of horizontal dispersion");
    itsAttr[DY] = Attributes::makeReal
                  ("DY", "Initial vertical dispersion");
    itsAttr[DPY] = Attributes::makeReal
                   ("DPY", "Initial slope of vertical dispersion");
    itsAttr[XC] = Attributes::makeReal
                  ("XC", "Initial horizontal position");
    itsAttr[PXC] = Attributes::makeReal
                   ("PXC", "Initial horizontal slope");
    itsAttr[YC] = Attributes::makeReal
                  ("YC", "Initial vertical position");
    itsAttr[PYC] = Attributes::makeReal
                   ("PYC", "Initial vertical slope");
    itsAttr[TC] = Attributes::makeReal
                  ("TC", "Initial time difference");
    itsAttr[PTC] = Attributes::makeReal
                   ("PTC", "Initial relative momentum error");

    // READ ONLY.
    itsAttr[LENGTH] = Attributes::makeReal
                      ("LENGTH", "Total length in m");
    itsAttr[LENGTH].setReadOnly(true);

    itsAttr[MU1] = Attributes::makeReal
                   ("MU1", "Phase for mode 1");
    itsAttr[MU1].setReadOnly(true);

    itsAttr[MU2] = Attributes::makeReal
                   ("MU2", "Phase for mode 2");
    itsAttr[MU2].setReadOnly(true);

    itsAttr[MU3] = Attributes::makeReal
                   ("MU3", "Phase for mode 3");
    itsAttr[MU3].setReadOnly(true);

    itsAttr[DELTAP] = Attributes::makeReal
                      ("DELTAP", "Differential momentum variation");
    itsAttr[DELTAP].setReadOnly(true);
}


Insertion::Insertion(const string &name, Insertion *parent):
    Twiss(name, parent)
{}


Insertion::~Insertion()
{}


Insertion *Insertion::clone(const string &name) {
    return new Insertion(name, this);
}


void Insertion::fill() {
    if(refill) {
        // Set truncation order.
        FTps<double, 6>::setGlobalTruncOrder(order);

        // Reset the mapper.
        if(itsAttr[INIT]) {
            TableRowRep rowrep = Attributes::getTableRow(itsAttr[INIT]);
            Table *init = Table::find(rowrep.getTabName());

            if(matches(init)) {
                PlaceRep pinit = rowrep.getPosition();
                Twiss   *tinit = dynamic_cast<Twiss *>(init);
                Row     &row   = tinit->findRow(pinit);
                tinit->fill();
                orbit = tinit->getOrbit(row);
                curly_A = tinit->getCurlyA(row);
                itsMapper->setMap(LinearMap<double, 6>() + orbit);
            } else {
                throw OpalException("Insertion::fill()",
                                    "Table \"" + rowrep.getTabName() +
                                    "\" is not suitable for initialising insertion \"" +
                                    getOpalName() + "\".");
            }
        } else {
            FMatrix<double, 6, 6> B, H;
            for(int i = 0; i < 6; ++i) B(i, i) = H(i, i) = 1.0;

            // Betatron matrix.
            B(0, 0) = sqrt(Attributes::getReal(itsAttr[BETX]));
            B(1, 1) = 1.0 / B(0, 0);
            B(1, 0) = - Attributes::getReal(itsAttr[ALFX]) * B(1, 1);
            B(2, 2) = sqrt(Attributes::getReal(itsAttr[BETY]));
            B(3, 3) = 1.0 / B(2, 2);
            B(3, 2) = - Attributes::getReal(itsAttr[ALFY]) * B(3, 3);
            B(4, 4) = 1.0;
            B(5, 5) = 1.0;

            // Dispersion matrix.
            H(0, 5) = Attributes::getReal(itsAttr[DX]);
            H(1, 5) = Attributes::getReal(itsAttr[DPX]);
            H(2, 5) = Attributes::getReal(itsAttr[DY]);
            H(3, 5) = Attributes::getReal(itsAttr[DPY]);
            H(4, 0) = - H(1, 5);
            H(4, 1) =   H(0, 5);
            H(4, 2) = - H(3, 5);
            H(4, 3) =   H(2, 5);

            curly_A = H * B;
            orbit(0) = Attributes::getReal(itsAttr[XC]);
            orbit(1) = Attributes::getReal(itsAttr[PXC]);
            orbit(2) = Attributes::getReal(itsAttr[YC]);
            orbit(3) = Attributes::getReal(itsAttr[PYC]);
            orbit(4) = Attributes::getReal(itsAttr[TC]);
            orbit(5) = Attributes::getReal(itsAttr[PTC]);
            itsMapper->setMap(LinearMap<double, 6>() + orbit);
        }

        put();
        refill = false;

        // Fill in the read-only data.
        const Row &row = itsTable->back();
        double arc = getS(row);
        Attributes::setReal(itsAttr[LENGTH], arc);
        Attributes::setReal(itsAttr[MU1], getMUi(row, 0));
        Attributes::setReal(itsAttr[MU2], getMUi(row, 1));
        Attributes::setReal(itsAttr[MU3], getMUi(row, 2));
    }
}


void Insertion::makeTFS(std::ostream &os, const CellArray &cells) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(12);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Write table descriptors.
    os << "@ TYPE     %s  TWISSTRACK\n";
    tfsTwissDescriptors(os);
    tfsSummary(os);

    // Write table body.
    tfsBody(os, cells);

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}


void Insertion::printTable(std::ostream &os, const CellArray &cells) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Print table header.
    printTableTitle(os, "Track lattice functions");

    // Print table body.
    printTableBody(os, cells);

    // Write table specific summary.
    const Row &row = itsTable->back();
    os << "Insert length =  " << setw(16) << getS(row)
       << "    Mux =        " << setw(16) << getMUi(row, 0)
       << "    Muy =        " << setw(16) << getMUi(row, 1)
       << '\n'
       << "                                 "
       << "    BetaX(max) = " << setw(16) << Attributes::getReal(itsAttr[BETXMAX])
       << "    BetaY(max) = " << setw(16) << Attributes::getReal(itsAttr[BETYMAX]) << '\n'
       << "                                 "
       << "    x(max) =     " << setw(16) << Attributes::getReal(itsAttr[XCMAX])
       << "    y(max) =     " << setw(16) << Attributes::getReal(itsAttr[YCMAX]) << '\n'
       << "                                 "
       << "    x(rms) =     " << setw(16) << Attributes::getReal(itsAttr[XCRMS])
       << "    y(rms) =     " << setw(16) << Attributes::getReal(itsAttr[YCRMS]) << '\n'
       << "                                 "
       << "    Dx(max) =    " << setw(16) << Attributes::getReal(itsAttr[DXMAX])
       << "    Dy(max) =    " << setw(16) << Attributes::getReal(itsAttr[DYMAX]) << '\n'
       << "                                 "
       << "    Dx(rms) =    " << setw(16) << Attributes::getReal(itsAttr[DXRMS])
       << "    Dy(rms) =    " << setw(16) << Attributes::getReal(itsAttr[DYRMS]) << '\n';

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}


void Insertion::tfsTwissDescriptors(std::ostream &os) const {
    tfsTableDescriptors(os);
    os << "@ LINE     %s  " << itsLine << '\n';
    beam->tfsDescriptors(os);
    os << "@ LENGTH   %le " << Attributes::getReal(itsAttr[LENGTH]) << '\n'
       << "@ MU1      %le " << Attributes::getReal(itsAttr[MU1]) << '\n'
       << "@ MU2      %le " << Attributes::getReal(itsAttr[MU2]) << '\n'
       << "@ MU3      %le " << Attributes::getReal(itsAttr[MU3]) << '\n';
}
