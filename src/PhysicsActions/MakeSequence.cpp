// ------------------------------------------------------------------------
// $CVSfile: MakeSequence.cc,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MakeSequence
//   The class for the OPAL MAKESEQ command.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:22:04 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "PhysicsActions/MakeSequence.h"
#include "AbsBeamline/Drift.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/ObjectFunction.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Algorithms/DefaultVisitor.h"
#include "Attributes/Attributes.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <fstream>

extern Inform *gmsg;

// Class MakeSequence
// ------------------------------------------------------------------------

namespace  MakeSequenceNS {

    // Functor for flagging an object.
    struct ObjectFlagger: ObjectFunction {
        virtual void operator()(Object *) const;
    };

    void ObjectFlagger::operator()(Object *object) const {
        // Only output objects which have a parent, and which are not built-in.
        object->setFlag(object->getParent() != 0  &&  ! object->isBuiltin());
    }

    // Functor for saving an element.
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
                *gmsg << object;
            }
            object->setFlag(false);
        }
    }

    // Functor for saving a variable.
    struct VariableWriter: ObjectFunction {
        VariableWriter(std::ostream &ostr): os(ostr) { }
        virtual void operator()(Object *) const;
    private:
        std::ostream &os;
    };

    void VariableWriter::operator()(Object *object) const {
        if(object->isFlagged() && dynamic_cast<ValueDefinition *>(object)) {
            *gmsg << object;
            object->setFlag(false);
        }
    }

    // Visitor class for writing the sequence.
    class SequenceWriter: public DefaultVisitor {
    public:
        // Construction/destruction.
        SequenceWriter(const Beamline &beamline, const std::string &name,
                       std::ostream &os);
        virtual ~SequenceWriter();

        // Override.
        virtual void execute();

        // Visit drift; must override to do nothing.
        virtual void visitDrift(const Drift &);

    protected:
        // Apply the default to an element.
        // All visitXXX() methods use applyDefault() which is overridden here.
        virtual void applyDefault(const ElementBase &element);

    private:
        // The sequence name.
        std::string itsName;

        // The output stream to be written.
        std::ostream &itsStream;

        // The accumulated length.
        double sum_length;
    };

    SequenceWriter::SequenceWriter(const Beamline &beamline, const std::string &name,
                                   std::ostream &os):
        DefaultVisitor(beamline, false, false),
        itsName(name), itsStream(os)
    {}

    void SequenceWriter::execute() {
        std::string comment = "// ";
        std::string line(72, '-');
        itsStream << comment << line << '\n'
                  << comment << "Sequence definition.\n"
                  << comment << line << '\n'
                  << itsName << ":SEQUENCE,REFER=CENTRE";
        itsStream << ",L=" << sum_length << ";\n";

        sum_length = 0.0;
        DefaultVisitor::execute();
        itsStream << "ENDSEQUENCE;" << std::endl;

    }

    SequenceWriter::~SequenceWriter()
    {}

    void SequenceWriter::visitDrift(const Drift &drift) {
        sum_length += drift.getElementLength();
    }

    void SequenceWriter::applyDefault(const ElementBase &element) {
        std::string objectName = element.getName();
        if(objectName[0] != '#') {
            Element *elem = Element::find(objectName);
            sum_length -= elem->getEntrance(Element::IS_CENTRE);
            itsStream << "   " << objectName << ",AT=" << sum_length;
            itsStream << ";\n";
            sum_length += elem->getExit(Element::IS_CENTRE);
        }
    }

    // The attributes of class MakeSequence.
    enum {
        LINE,        // The lattice to be used.
        NAME,        // The name for the new sequence.
        FNAME,       // The file to be written.
        SIZE
    };
}

using namespace MakeSequenceNS;

MakeSequence::MakeSequence():
    Action(SIZE, "MAKESEQ",
           "The \"MAKESEQ\" statement constructs a flat sequence from a "
           "\"LINE\" object.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of the lattice to be flattened");
    itsAttr[NAME] = Attributes::makeString
                    ("NAME",
                     "Name to be given to the generated seqence (default = original name).");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE",
                      "Name to be given to the generated file (default = new sequence name).");
}


MakeSequence::MakeSequence(const std::string &name, MakeSequence *parent):
    Action(name, parent)
{}


MakeSequence::~MakeSequence()
{}


MakeSequence *MakeSequence::clone(const std::string &name) {
    return new MakeSequence(name, this);
}


void MakeSequence::execute() {
    // Get relevant names.
    std::string useName = Attributes::getString(itsAttr[LINE]);
    std::string name = Attributes::getString(itsAttr[NAME]);
    std::string file = Attributes::getString(itsAttr[FNAME]);
    if(name.empty()) name = useName;
    if(file.empty()) file = name;

    // Find BeamSequence definition.
    BeamSequence *use = BeamSequence::find(useName);

    // Open the output stream.
    std::ofstream os(file.c_str());
    if(os.bad()) {
        throw OpalException("MakeSequence::execute()",
                            "Unable to open output stream \"" + file + "\".");
    }
    os.precision(12);
    std::string line(72, '-');
    std::string comment = "// ";

    // Flag all objects which should be output.
    OpalData::getInstance()->apply(MakeSequenceNS::ObjectFlagger());

    // Write all variables.
    os << comment << line << '\n'
       << comment << "Variable definitions." << '\n'
       << comment << line << '\n';
    OpalData::getInstance()->apply(MakeSequenceNS::VariableWriter(os));
    os << '\n';

    // Write all elements.
    os << comment << line << '\n'
       << comment << "Element definitions." << '\n'
       << comment << line << '\n';
    OpalData::getInstance()->apply(MakeSequenceNS::ElementWriter(os));
    os << std::endl;

    // Write the sequence.

    SequenceWriter writer(*use->fetchLine(), name, os);
    writer.execute();
    os.close();
}
