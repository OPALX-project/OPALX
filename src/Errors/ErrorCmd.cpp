// ------------------------------------------------------------------------
// $RCSfile: ErrorCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorCmd
//   The class for the OPAL ERROR command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorCmd.h"
#include "AbstractObjects/BeamSequence.h"
#include "Algorithms/Flagger.h"
#include "Attributes/Attributes.h"
#include "Errors/AlignHandler.h"
#include "Errors/AlignRemover.h"
#include "Errors/Error.h"
#include "Errors/ErrorParser.h"
#include "Errors/MPHandler.h"
#include "Errors/MPRemover.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"


// Class ErrorCmd
// ------------------------------------------------------------------------

// The attributes of class ErrorCmd.
namespace {
  enum {
    LINE,      // The lattice to be assigned field errors.
    SEED,      // The new error seed, if given.
    ADD,       // If true, errors are additive.
    CLEAR,     // If true, clear selection flags first.
    SIZE
  };
}


ErrorCmd::ErrorCmd():
  Action(SIZE, "ERROR",
	 "The \"ERROR\" command initiates error definition mode.")
{
  itsAttr[LINE] = Attributes::makeString
    ("LINE", "Beam line or sequence for which errors must be generated");
  itsAttr[SEED] = Attributes::makeReal
    ("SEED", "If entered, sets the seed for the random generator");
  itsAttr[ADD] = Attributes::makeBool
    ("ADD", "If true, error definitions are added on top of existing ones");
  itsAttr[CLEAR] = Attributes::makeBool
    ("CLEAR", "If true, all selections are cleared before execution.");
}


ErrorCmd::ErrorCmd(const string &name, ErrorCmd *parent):
  Action(name, parent)
{}


ErrorCmd::~ErrorCmd()
{}


ErrorCmd *ErrorCmd::clone(const string &name)
{
  return new ErrorCmd(name, this);
}


void ErrorCmd::execute()
{
  Error::block = new Error();

  // Define the active line.
  const string &lineName = Attributes::getString(itsAttr[LINE]);
  Error::block->itsLine = BeamSequence::find(lineName)->fetchLine();

  // If "CLEAR", remove all errors.
  if (Attributes::getBool(itsAttr[CLEAR])) {
    {
      AlignRemover remover;
      AlignHandler visitor(*Error::block->itsLine, remover, true);
      visitor.execute();
    }
    {
      MPRemover remover;
      MPHandler visitor(*Error::block->itsLine, remover, true);
      visitor.execute();
    }
  }

  // Deselect all elements.
  Flagger flagger(*Error::block->itsLine, false);
  flagger.execute();

  // Store the add flag.
  Error::block->addError = Attributes::getBool(itsAttr[ADD]);

  // If given, set the random seed.
  if (itsAttr[SEED]) {
    Options::seed = int(Round(Attributes::getReal(itsAttr[SEED])));
    Options::rangen.init55(Options::seed);
  }

  // Run the error parser.
  Error::block->parser.run();

  // Deselect all elements (re-use the flagger declared above).
  flagger.execute();

  // Clean up.
  delete Error::block;
  Error::block = 0;
}
