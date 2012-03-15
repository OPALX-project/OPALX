// ------------------------------------------------------------------------
// $RCSfile: SetIntegrator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SetIntegrator
//   The class for the OPAL SETINTEGRATOR command.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/03/28 21:27:54 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "PhysicsActions/SetIntegrator.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/Integrator.h"
#include "AbsBeamline/Multipole.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/MPSplitIntegrator.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Utilities/Options.h"
#include <iostream>

using std::cerr;


// Class SetIntegrator
// ------------------------------------------------------------------------

namespace {

    // Visitor class Setter for attaching the integrators.
    class Setter: public DefaultVisitor {

    public:
        // Construction/destruction.
        Setter(const Beamline &, const string &type, int slices);
        virtual ~Setter();

        // Override method for integratorconst  visit.
        virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

        // Return error counts.
        void getCounts(int &add, int &rem, int &rep) const;

    private:
        // The type of integrator.
        const string itsType;

        // The number of slices.
        int itsSlices;

        // Keeps track of inserted integrators.
        int added, removed, replaced;
    };


    Setter::Setter(const Beamline &beamline, const string &type, int slices):
        DefaultVisitor(beamline, false, false),
        itsType(type), itsSlices(slices), added(0), removed(0), replaced(0)
    {}


    Setter::~Setter()
    {}


    void Setter::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
        // The pointer "base" points to the actual element or integrator.
        ElementBase *temp = fep.getElement();
        ElementBase *base = temp->removeAlignWrapper();

        // The pointer "wrap" points to an AlignWrapper, if it exist,
        // otherwise it is NULL.
        AlignWrapper *wrap = dynamic_cast<AlignWrapper *>(temp);

        if(dynamic_cast<Beamline *>(base)) {
            // Enter beam line via the default action.
            DefaultVisitor::visitFlaggedElmPtr(fep);
        } else if(fep.getSelectionFlag()) {
            // Single element is selected.
            if(Integrator *i = dynamic_cast<Integrator *>(base)) {
                // The embedded element should be a multipole.
                if(Multipole *mult = dynamic_cast<Multipole *>(i->getElement())) {
                    // Replace/remove an existing integrator.
                    ElementBase *mpsi;
                    if(itsType == "NONE") {
                        // Going to link the multipole directly, without an integrator.
                        removed++;
                        mpsi = mult;
                    } else {
                        // Going to link the multipole with a new integrator.
                        replaced++;
                        mpsi = new MPSplitIntegrator(mult, itsSlices);
                    }
                    if(wrap) {
                        wrap->setElement(mpsi);
                    } else {
                        const_cast<FlaggedElmPtr &>(fep).setElement(mpsi);
                    }
                }
            } else if(Multipole *mult = dynamic_cast<Multipole *>(base)) {
                // Going to wrap the multipole with a new integrator.
                ElementBase *mpsi = new MPSplitIntegrator(mult, itsSlices);
                added++;

                if(wrap) {
                    wrap->setElement(mpsi);
                } else {
                    const_cast<FlaggedElmPtr &>(fep).setElement(mpsi);
                }
            }
            // Nothing to be done for any other element.
        }
    }


    void Setter::getCounts(int &add, int &rem, int &rep) const {
        add = added;
        rem = removed;
        rep = replaced;
    }


    // The attributes of class SetIntegrator.
    enum {
        LINE,        // The lattice to be used.
        TYPE,        // The type of integrator.
        SLICES,      // The number of slices.
        SIZE
    };
}


SetIntegrator::SetIntegrator():
    Action(SIZE, "SETINT",
           "The \"SETINT\" statement attaches special integrators to selected "
           "elements.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of the lattice to be changed");
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Type of integrator to be attached: default = THIN", "THIN");
    itsAttr[SLICES] = Attributes::makeReal
                      ("SLICES", "Number of slices to be used for \"TYPE=THIN\"", 1.0);
}


SetIntegrator::SetIntegrator(const string &name, SetIntegrator *parent):
    Action(name, parent)
{}


SetIntegrator::~SetIntegrator()
{}


SetIntegrator *SetIntegrator::clone(const string &name) {
    return new SetIntegrator(name, this);
}


void SetIntegrator::execute() {
    // Find BeamSequence definition.
    BeamSequence *use = BeamSequence::find(Attributes::getString(itsAttr[LINE]));

    // Initialize the visitor.
    int slices = int(Attributes::getReal(itsAttr[SLICES]));
    Setter splitter(*use->fetchLine(), Attributes::getString(itsAttr[TYPE]), slices);
    splitter.execute();

    // Confirm action to user.
    if(Options::info) {
        int add = 0, rem = 0, rep = 0;
        splitter.getCounts(add, rem, rep);
        cerr << '\n';

        if(rem == 0) {
            cerr << "No integrator removed.\n";
        } else if(rem == 1) {
            cerr << "1 integrator removed.\n";
        } else if(rem > 1) {
            cerr << rem << " integrators removed.\n";
        }

        if(rep == 0) {
            cerr << "No integrator replaced.\n";
        } else if(rep == 1) {
            cerr << "1 integrator replaced.\n";
        } else if(rep > 1) {
            cerr << rep << " integrators replaced.\n";
        }

        if(add == 0) {
            cerr << "No integrator installed.\n";
        } else if(add == 1) {
            cerr << "1 integrator installed.\n";
        } else if(add > 1) {
            cerr << add << " integrators installed.\n";
        }

        cerr << std::endl;
    }
}
