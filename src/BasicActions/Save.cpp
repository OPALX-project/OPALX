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
#include "config.h"

extern Inform *gmsg;

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
                *gmsg << object;//->print(*gmsg);
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
            *gmsg << object;//->print(*gmsg);
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
            *gmsg << object;//->print(*gmsg);
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


Save::Save(const std::string &name, Save *parent):
    Action(name, parent)
{}


Save::~Save()
{}


Save *Save::clone(const std::string &name) {
    return new Save(name, this);
}


void Save::execute() {
    std::string file = Attributes::getString(itsAttr[0]);
    std::ofstream os(file.c_str());

    if(os.bad()) {
        throw OpalException("Save::execute()",
                            "Unable to open output stream \"" + file + "\".");
    } else {
        // Flag all objects to be saved.
        OpalData::getInstance()->apply(ObjectFlagger());


        // Now save all objects according to categories.
        //JMJ adding some comment tags to saved output 25/10/2000
        //JMJ more of those 18/12/2000

        std::string comchar = "// ";

        os << comchar << "<OPAL Version " << PACKAGE_VERSION << " SVN version "
           << SVN_VERSION << "  (c) PSI, http://amas.web.psi.ch"
           << std::endl << ";" << std::endl ;

        os << comchar << "<Parameter definitions> ;" << std::endl ;
        OpalData::getInstance()->apply(ParameterWriter(os));
        os << comchar << "</Parameter definitions> ;"
           << std::endl << ";" << std::endl ;

        os << comchar << "<Element definitions> ;" << std::endl ;
        OpalData::getInstance()->apply(ElementWriter(os));
        os << comchar << "</Element definitions> ;"
           << std::endl << ";" << std::endl ;

        os << comchar << "<Line (and split element) definitions> ;"
           << std::endl ;

        os << comchar << "</Line (and split element) definitions> ;"
           << std::endl << ";" << std::endl ;

        os << comchar << "<Special definitions> ;" << std::endl ;
        OpalData::getInstance()->apply(SpecialWriter(os));
        os << comchar << "</Special definitions> ;"
           << std::endl << ";" << std::endl ;

        os << comchar << "<OPAL Version " << PACKAGE_VERSION << " SVN version "
           << SVN_VERSION << "  (c) PSI, http://amas.web.psi.ch"
           << std::endl << ";" << std::endl ;
    }
}


void Save::parse(Statement &statement) {
    parseShortcut(statement);
}
