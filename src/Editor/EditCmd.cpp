// ------------------------------------------------------------------------
// $RCSfile: EditCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditCmd
//   The class for the OPAL sequence editor SEQEDIT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditCmd.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Editor/Edit.h"
#include "Utilities/OpalException.h"


// Class EditCmd
// ------------------------------------------------------------------------

EditCmd::EditCmd():
    Action(1, "SEQEDIT",
           "The \"SEQEDIT\" command specifies a sequence for editing.") {
    itsAttr[0] = Attributes::makeString
                 ("SEQUENCE", "Name of sequence to be edited");
}


EditCmd::EditCmd(const std::string &name, EditCmd *parent):
    Action(name, parent)
{}


EditCmd::~EditCmd()
{}


EditCmd *EditCmd::clone(const std::string &name) {
    return new EditCmd(name, this);
}


void EditCmd::execute() {
    // The parser strategy for sequence edit mode.
    std::string seqName = Attributes::getString(itsAttr[0]);

    if(seqName == "") {
        throw OpalException("EditCmd::execute()",
                            "Sequence name is required for \"SEQEDIT\".");
    } else {
        // find sequence definition
        Object *object = OpalData::getInstance()->find(seqName);
        Sequence *seq = dynamic_cast<Sequence *>(object);

        if(seq == 0) {
            throw OpalException("EditCmd::execute()",
                                "Sequence \"" + seqName + "\" not found.");
        } else {
            // Execute edit block.
            Edit::block = new Edit(seq);
            Edit::block->parser.run();

            // Clean up.
            delete Edit::block;
            Edit::block = 0;
        }
    }
}
