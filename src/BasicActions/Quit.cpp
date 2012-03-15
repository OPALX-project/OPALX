// ------------------------------------------------------------------------
// $RCSfile: Quit.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Quit
//   The class for the OPAL QUIT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Quit.h"


// Class Quit
// ------------------------------------------------------------------------

Quit::Quit(): Action(0, "QUIT",
		     "The \"QUIT\" statement terminates program execution "
		     "or reading of a called file.")
{}


Quit::Quit(const string &name, Quit *parent):
  Action(name, parent)
{}


Quit::~Quit()
{}


Quit *Quit::clone(const string &name)
{
  return new Quit(name, this);
}


void Quit::execute()
{}
