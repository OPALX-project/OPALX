// ------------------------------------------------------------------------
// $RCSfile: Save.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Save
//   The base class for the OPAL SAVE command.
//
//   Note by JMJ 6/4/2000:
//   According to some old version of the manual this should have a PATTERNS option
//   to allow only a certain set of parameters to be saved.
//   Not so, it seems to just save everything.   I opale changes to the print method
//   of the ConcreteVar class to partly compensate for this deficiency: that allows
//   you to get the variables used in matching in OPAL input syntax, at least on the
//   main output file.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/14 07:02:44 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "BasicActions/Save.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Definition.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/ObjectFunction.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <fstream>
#include <string>


// Functors for flagging objects and saving special categories.
// ------------------------------------------------------------------------

namespace  SaveNS {

    // Functor for flagging an object.
    // ----------------------------------------------------------------------
    struct ObjectFlagger: ObjectFunction {
        virtual void operator()(Object *) const;
    };

    void ObjectFlagger::operator()(Object *object) const {
        // Only output objects which have a parent, and which are not built-in.
        object->setFlag(object->getParent() != 0  &&  ! object->isBuiltin());
    }

    // Functor for saving an element.
    // ----------------------------------------------------------------------
    struct ElementWriter: ObjectFunction {
        ElementWriter(std::ostream &ostr): os(ostr) { }
        virtual void operator()(Object *) const;
    private:
        std::ostream &os;
    };

    void ElementWriter::operator()(Object *object) const {
        if(object->isFlagged() && dynamic_cast<Element *>(object) &&
           ! dynamic_cast<BeamSequence *>(object)) {
            if(object->getOpalName()[0] != '#') {
                (*this)(object->getParent());
                object->print();
                object->setFlag(false);
            }
        }
    }

    // Functor for saving a parameter.
    // ----------------------------------------------------------------------
    struct ParameterWriter: ObjectFunction {
        ParameterWriter(std::ostream &ostr): os(ostr) { }
        virtual void operator()(Object *) const;
    private:
        std::ostream &os;
    };

    void ParameterWriter::operator()(Object *object) const {
        if(object->isFlagged() && dynamic_cast<ValueDefinition *>(object)) {
            object->print();
            object->setFlag(false);
        }
    }

    // Functor for saving a line or a sequence.
    // ----------------------------------------------------------------------
    struct LineWriter: ObjectFunction {
        LineWriter(std::ostream &ostr): os(ostr) { }
        virtual void operator()(Object *) const;
    private:
        std::ostream &os;
    };

    void LineWriter::operator()(Object *object) const {
        if(object->isFlagged() && dynamic_cast<BeamSequence *>(object)) {
            (*this)(object->getParent());
            object->print();
            object->setFlag(false);
        }
    }

    // Functor for saving a special definition.
    // ----------------------------------------------------------------------
    struct SpecialWriter: ObjectFunction {
        SpecialWriter(std::ostream &ostr): os(ostr) { }
        virtual void operator()(Object *) const;
    private:
        std::ostream &os;
    };

    void SpecialWriter::operator()(Object *object) const {
        if(object->isFlagged() && dynamic_cast<Definition *>(object)) {
            (*this)(object->getParent());
            object->print();
            object->setFlag(false);
        }
    }
}


using namespace SaveNS;



// Class Save
// ------------------------------------------------------------------------

Save::Save():
    Action(1, "SAVE",
           "The \"SAVE\" statement prints a list of all definitions,\n"
           "starting with constants, variables, and vectors,"
           "followed by elements, and finally all sequences.") {
    itsAttr[0] = Attributes::makeString
                 ("FILE", "Name of file to be written", "SAVE");
}


Save::Save(const string &name, Save *parent):
    Action(name, parent)
{}


Save::~Save()
{}


Save *Save::clone(const string &name) {
    return new Save(name, this);
}


void Save::execute() {
    string file = Attributes::getString(itsAttr[0]);
    std::ofstream os(file.c_str());

    if(os.bad()) {
        throw OpalException("Save::execute()",
                            "Unable to open output stream \"" + file + "\".");
    } else {
        // Flag all objects to be saved.
        OPAL.apply(ObjectFlagger());


        // Now save all objects according to categories.
        //JMJ adding some comment tags to saved output 25/10/2000
        //JMJ more of those 18/12/2000

        string comchar = "// ";
        int version = 9;
        if(Options::mad8)  {comchar = "!! " ; version = 8;}

        os << comchar << "<OPAL Version " << version
           << " data, created by OPAL Version 9.> ;"
           << endl << ";" << endl ;

        os << comchar << "<Parameter definitions> ;" << endl ;
        OPAL.apply(ParameterWriter(os));
        os << comchar << "</Parameter definitions> ;"
           << endl << ";" << endl ;

        os << comchar << "<Element definitions> ;" << endl ;
        OPAL.apply(ElementWriter(os));
        os << comchar << "</Element definitions> ;"
           << endl << ";" << endl ;

        os << comchar << "<Line (and split element) definitions> ;"
           << endl ;
        OPAL.apply(LineWriter(os));
        os << comchar << "</Line (and split element) definitions> ;"
           << endl << ";" << endl ;

        os << comchar << "<Special definitions> ;" << endl ;
        OPAL.apply(SpecialWriter(os));
        os << comchar << "</Special definitions> ;"
           << endl << ";" << endl ;

        os << comchar << "</OPAL Version " << version
           << " data, created by OPAL Version 9.> ;"
           << endl << ";" << endl ;

    }
}


void Save::parse(Statement &statement) {
    parseShortcut(statement);
}
