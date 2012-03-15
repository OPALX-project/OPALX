#ifndef OPAL_ErrorSave_HH
#define OPAL_ErrorSave_HH

// ------------------------------------------------------------------------
// $RCSfile: ErrorSave.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorSave
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Errors/AlignBase.h"
#include "Errors/MPBase.h"
#include <fstream>

class AlignWrapper;


// Class ErrorSave
// ------------------------------------------------------------------------
/// The ESAVE command.

class ErrorSave:
  public Action, private AlignBase, private MPBase {

public:

  /// Exemplar constructor.
  ErrorSave();

  virtual ~ErrorSave();

  /// Make clone.
  virtual ErrorSave *clone(const string &name);

  /// Execute the command.
  virtual void execute();

  /// Save error field. 
  //  Use the name and occurrence count for the current element.  Relative
  //  errors are related to [b]mainField[/b]. The field error is fetched
  //  from [b]errorField[/b].  This class's execute() method uses a MPHandler
  //  to call this method for all elements which have a field error.
  virtual void fieldError(const string &name, int occur,
			  const BMultipoleField &mainField,
			  BMultipoleField &errorField);

  /// Save misalignement.
  //  All information is taken from the given AlignWrapper.
  //  This class's execute() method uses an AlignHandler to call this
  //  method for all elements which have a misalignment.
  virtual void misalignment(const AlignWrapper &, int occur);

private:

  // Not implemented.
  ErrorSave(const ErrorSave &);
  void operator=(const ErrorSave &);

  // Clone constructor.
  ErrorSave(const string &name, ErrorSave *parent);

  // Output stream; "mutable" is required to enable writing.
  mutable std::ofstream os;

  // Data for command execution.
  bool doAlign, doField;
  double itsRadius;
  int itsOrder;
};

#endif // OPAL_ErrorSave_HH
