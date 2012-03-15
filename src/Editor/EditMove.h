#ifndef OPAL_EditMove_HH
#define OPAL_EditMove_HH

// ------------------------------------------------------------------------
// $RCSfile: EditMove.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EditMove
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Editor.h"


// Class EditMove
// ------------------------------------------------------------------------
/// The sequence editor MOVE command.

class EditMove: public Editor {

public:

  /// Exemplar constructor.
  EditMove();

  virtual ~EditMove();

  /// Make clone.
  virtual EditMove *clone(const string &name);

  /// Execute the command.
  virtual void execute();

private:

  // Not implemented.
  EditMove(const EditMove &);
  void operator=(const EditMove &);

  // Clone constructor.
  EditMove(const string &name, EditMove *parent);
};

#endif // OPAL_EditMove_HH
