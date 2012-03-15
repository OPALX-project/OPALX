#ifndef OPAL_Flatten_HH
#define OPAL_Flatten_HH

// ------------------------------------------------------------------------
// $RCSfile: Flatten.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Flatten<T>
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Tables/RangeSelector.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "BeamlineCore/PatchRep.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Beamlines/TBeamline.h"
#include "Beamlines/FlaggedElmPtr.h"

class RangeRep;


// Template class Flatten
// ------------------------------------------------------------------------
/// Flatten a beamline.
//  The type of beam line members in the flat line is given as a template
//  parameter.  It may be any class derived from ElmPtr.  Execution of
//  this visitor yields a flat line of the type defined by the template
//  parameter Member.

template <class Member> class Flatten: public RangeSelector {

public:

  /// Constructor.
  //  Attach this visitor to [b]bl[/b], [b]mem[/b] will receive the range
  //  [b]range[/b] ot the flat line.
  Flatten(const Beamline &bl, TBeamline<Member> &mem, const RangeRep &range);

  virtual ~Flatten();

  /// Apply the algorithm to the top-level beamline.
  virtual void execute();

protected:

  /// The operation to be done for beamlines.
  virtual void handleBeamline(const FlaggedElmPtr &);

  /// The operation to be done for elements.
  virtual void handleElement(const FlaggedElmPtr &);

  /// The flat list to be filled.
  TBeamline<Member> &itsTable;

private:

  // Not implemented.
  Flatten();
  Flatten(const Flatten<Member> &);
  void operator=(const Flatten<Member> &);
};


// Implementation of template class Flatten
// ------------------------------------------------------------------------

template <class Member> inline
Flatten<Member>::Flatten(const Beamline &beamline, TBeamline<Member> &list,
			 const RangeRep &range):
  RangeSelector(beamline, range), itsTable(list)
{}


template <class Member> inline
Flatten<Member>::~Flatten()
{}


template <class Member> inline
void Flatten<Member>::execute()
{
  itsTable.push_back(Member(Element::find("#S")->getElement(), 0));
  RangeSelector::execute();
  itsTable.push_back(Member(Element::find("#E")->getElement(), 0));
}


template <class Member> inline void
Flatten<Member>::handleBeamline(const FlaggedElmPtr &fep)
{
  AlignWrapper &wrap =
    dynamic_cast<AlignWrapper &>(*fep.getElement()->makeAlignWrapper());

  if (itsRange.isActive()) {
    // Append the entrance misalignment patch.
    Euclid3D transform = wrap.getEntranceTransform();
    PatchRep *patch = new PatchRep("[BEGIN]" + wrap.getName());
    patch->setPatch(transform);
    Member member(patch, 0);
    itsTable.push_back(member);
  }

  // Delegate algorithm to the beamline.
  RangeSelector::handleBeamline(fep);

  if (itsRange.isActive()) {
    // Append the exit misalignment patch.
    Euclid3D transform = wrap.getExitTransform();
    PatchRep *patch = new PatchRep("[END]" + wrap.getName());
    patch->setPatch(transform);
    Member member(patch, 0);
    itsTable.push_back(member);
  }
}


template <class Member> void
Flatten<Member>::handleElement(const FlaggedElmPtr &fep)
{
  if (itsRange.isActive()) {
    const string &name = fep.getElement()->getName();
    if (name[0] != '#') {
      Member member(fep);
      itsTable.push_back(member);
    }
  }
}

#endif // OPAL_Flatten_HH
