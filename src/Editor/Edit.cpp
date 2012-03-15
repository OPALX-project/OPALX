// ------------------------------------------------------------------------
// $RCSfile: Edit.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Edit
//   This class is the work horse for OPAL sequence editor commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/12/15 10:13:47 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Editor/Edit.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Marker.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/RangeRep.h"
#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/Flagger.h"
#include "BeamlineCore/DriftRep.h"
#include "Beamlines/TBeamline.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Elements/OpalDrift.h"
#include "Lines/Sequence.h"
#include "Lines/SequenceMember.h"
#include "MemoryManagement/Pointer.h"
#include "Tables/Selector.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/RegularExpression.h"
#include <iostream>


namespace {

  // Class EditFlat
  //   Flatten the sequence being edited.
  // ----------------------------------------------------------------------

  class EditFlat: public DefaultVisitor {
    
  public:
    
    // Construction/destruction
    EditFlat(const Beamline &, Edit::TLine &, Sequence::ReferenceType);
    virtual ~EditFlat();

    
    //: Apply the algorithm to a drift.
    //  Accumulate drift length and skip.
    virtual void visitDrift(const Drift &);
    
    //: Apply the algorithm to a marker.
    //  Suppress begin and end markers.
    virtual void visitMarker(const Marker &);
    
  private:
    
    // Not implemented.
    EditFlat();
    EditFlat(const EditFlat &);
    void operator=(const EditFlat &);

    // Apply default:
    // Append element to flat sequence.
    virtual void applyDefault(const ElementBase &);
    
    // The line to be filled during application of this visitor.
    Edit::TLine &newLine;

    // The accumulated total length.
    double sumLength;
    
    // The reference type for the outer sequence.
    Sequence::ReferenceType reference;
  };
  
  
  // Implementation of class EditFlat
  // ----------------------------------------------------------------------
  
  EditFlat::EditFlat(const Beamline &line, Edit::TLine &toLine,
		     Sequence::ReferenceType refer):
    DefaultVisitor(line, false, false), newLine(toLine), reference(refer)
    {
      sumLength = 0.0;
    }
  
  
  EditFlat::~EditFlat()
    {}
  
  
  void EditFlat::visitDrift(const Drift &drift)
    {
      // Accumulate drift length.
      sumLength += drift.getElementLength();
    }
  
  
  void EditFlat::visitMarker(const Marker &marker)
    {
      const string &name = marker.getName();
      if (name[0] != '#') applyDefault(marker);
    }
  
  
  void EditFlat::applyDefault(const ElementBase &element)
    {
      double length = element.getElementLength();
      SequenceMember member;
      Element *elem = Element::find(element.getName());
      member.setElement(elem->getElement());

      // ada 4.5 2000 to speed up matching, add a pointer to
      // opal elements in order to avoid serching the opal elements
      member.OpalElement = elem;

      switch (reference) {
      case Sequence::IS_ENTRY:
	member.itsPosition = sumLength;
	break;
      case Sequence::IS_CENTRE:
	member.itsPosition = sumLength + length / 2.0;
	break;
      case Sequence::IS_EXIT:
	member.itsPosition = sumLength + length;
	break;
      }
      
      sumLength += length;
      newLine.push_back(member);
    }
};


// Class Edit
// This code assumes that the sequence has at least 3 members,
// namely the begin and end markers, and a drift between the two.
// We only look at the user-defined positions, and skip the drifts.
// ------------------------------------------------------------------------

Edit *Edit::block = 0;


Edit::Edit(Sequence *seq):
  itsSequence(seq), itsLine(), isModified(false), parser()
{
  Pointer<Sequence> clone = seq->copy(seq->getOpalName());
  itsLine = clone->fetchLine();
  selectClear();
}


Edit::~Edit()
{}


bool Edit::cycle(const PlaceRep &init)
{
  // NOTE: We can cycle the top level only.
  // Get length of sequence.
  PlaceRep pos(init);
  double length = itsSequence->getLength();

  // Point "first" to first element (may be end marker).
  iterator first = itsLine->begin();
  ++first; ++first;

  // Point "last" to end marker.
  iterator last = itsLine->end();
  --last;

  for (iterator start = first; start != last; ) {
    pos.enter(*start);
    
    if (pos.isActive()) {
      // Update positions preceding "start".
      for (iterator i = first; i != start; ) {
	(*i).itsPosition += length;
	++i; ++i;
      }

      // Move all elements from "first" (included) to "start" (excluded),
      // and the drifts following each element to just before "last".
      itsLine->splice(last, *itsLine, first, start);
      
      // Point "first" to new start element.
      iterator first = itsLine->begin();
      ++first; ++first;

      // Point "last" to end marker.
      iterator last = itsLine->end();
      --last;

      // Update the lengths to start with zero.
      double startPos = (*first).itsPosition;
      for (iterator i = first; i != last; ) {
	(*i).itsPosition -= startPos;
	++i; ++i;
      }
      
      isModified = true;
      return true;
    }
    
    pos.leave(*start);
    start++; start++;
  }

  return false;
}


void Edit::finish(const string &newName)
{
  selectClear();

  if (newName == itsSequence->getOpalName()) {
    // If name is unchanged and sequence is modified,
    // replace the sequence line.
    if (isModified) {
      itsSequence->storeLine(*itsLine);
    }
  } else {
    // If name is changed, build new modified sequence.
    Sequence *model = dynamic_cast<Sequence *>(OPAL.find("SEQUENCE"));
    Pointer<Sequence> copy = model->clone(newName);
    copy->copyAttributes(*itsSequence);
    copy->storeLine(*itsLine);
    OPAL.define(&*copy);
    copy->update();
  }
  
  OPAL.makeDirty(&*itsSequence);
}


void Edit::flatten()
{
  Sequence::ReferenceType refer = itsSequence->getReference();
  TLine newLine;
  EditFlat flatten(*itsLine, newLine, refer);
  flatten.execute();
  itsSequence->addEndMarkers(newLine);
  itsSequence->insertDrifts(newLine);
  itsLine->swap(newLine);
  isModified = true;
}


int Edit::installMultiple(ElementBase *elm, double at)
{
  return installMultiple(false, *itsLine, elm, at);
}


int Edit::installSingle(const PlaceRep &pos, ElementBase *elm, double at)
{
  PlaceRep from(pos);
  return installSingle(false, *itsLine, from, elm, at);
}


int Edit::moveMultiple(double by)
{
  return moveMultiple(false, *itsLine, by);
}


int Edit::moveSingleAbs(const PlaceRep &pos, double to)
{
  PlaceRep it(pos);
  return moveSingleAbs(false, *itsLine, it, to);
}


int Edit::moveSingleRel(const PlaceRep &pos, const PlaceRep &from, double to)
{
  PlaceRep it(pos);
  PlaceRep frm(from);
  return moveSingleRel(false, *itsLine, it, frm, to);
}


void Edit::reflect()
{
  itsLine = reflect(*itsLine);
}


int Edit::removeMultiple()
{
  return removeMultiple(false, *itsLine);
}


int Edit::removeSingle(const PlaceRep &pos)
{
  PlaceRep it(pos);
  return removeSingle(false, *itsLine, it);
}


int Edit::replaceMultiple(ElementBase *elm)
{
  return replaceMultiple(false, *itsLine, elm);
}


int Edit::replaceSingle(const PlaceRep &pos, ElementBase *elm)
{
  PlaceRep it(pos);
  return replaceSingle(false, *itsLine, it, elm);
}


int Edit::select(const RangeRep &rng, const string &cls,
		 const string &typ, const string &patt)
{
  Selector sel(*itsLine, rng, cls, typ, patt);
  sel.execute();
  return sel.getCount();
}


void Edit::selectClear()
{
  Flagger flagger(*itsLine, false);
  flagger.execute();
}


void Edit::selectFull()
{
  Flagger flagger(*itsLine, true);
  flagger.execute();
}


void Edit::install(TLine &newList, ElementBase *elm, double at)
{
  SequenceMember member1;

  // Push the installed element.
  member1.setElement(elm->copyStructure());

  // ada 4.5 2000 to speed up matching, add a pointer to
  // opal elements in order to avoid serching the opal elements
  member1.OpalElement = Element::find(elm->getName());

  member1.itsType = SequenceMember::GLOBAL;
  member1.itsPosition = at;
  newList.push_back(member1);

  // Push a generated drift.
  SequenceMember member2;
  DriftRep *drift = new DriftRep();
  drift->setName("[DRIFT]");
  member2.setElement(drift);
  member2.itsType = SequenceMember::GENERATED;
  member2.itsPosition = 0.0;
  newList.push_back(member2);
}


int Edit::installMultiple(bool shared, TLine &bl, ElementBase *elm, double at)
{
  int count = 0;
  TLine newList;

  // Point "first" at begin marker.
  iterator first = bl.begin();

  // Point "last" at end marker.
  iterator last = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    if (i->getSelectionFlag()) {
      if (shared && Options::warn) {
	invalidShare("INSTALL");
      } else {
	install(newList, elm, i->itsPosition + at);
	++count;
      }
    } else {
      ElementBase *base = i->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("INSTALL");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += installMultiple(shared || base->isSharable(), *line, elm, at);
      }
    }

    ++i; ++i;
  }

  merge(bl, newList);
  return count;
}


int Edit::installSingle
(bool shared, TLine &bl, PlaceRep &it, ElementBase *elm, double at)
{
  int count = 0;
  TLine newList;

  // Point "first" at begin marker.
  iterator first = bl.begin();

  // Point "last" at end marker.
  iterator last = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    it.enter(*i);
    
    if (it.isActive()) {
      if (shared && Options::warn) {
	invalidShare("INSTALL");
      } else {
	install(newList, elm, i->itsPosition + at);
	++count;
      }
    } else {
      ElementBase *base = i->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("INSTALL");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += installSingle(shared || base->isSharable(),
			       *line, it, elm, at);
      }
    }
    
    it.leave(*i);
    ++i; ++i;
  }

  merge(bl, newList);
  return count;
}


void Edit::merge(TLine &bl, TLine &newList)
{
  // Point "first" at first element.
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last at end marker.
  iterator last = bl.end();
  --last;

  // For all positions in original line.
  for (iterator i = first; i != last; ) {
    double sOld = i->itsPosition;

    // Insert new positions which belong before the current position.
    iterator j = newList.begin();
    iterator k = j;
    while (k != newList.end()  &&  k->itsPosition < sOld) {
      ++k; ++k;
    }

    if (j != k) {
      bl.splice(i, newList, j, k);
      isModified = true;
    }

    ++i; ++i;
  }

  // Insert any remaining positions before end marker.
  // This happens even if their position is beyond the end marker.
  if (! newList.empty()) {
    bl.splice(last, newList, newList.begin(), newList.end());
    isModified = true;
  }
}


int Edit::moveMultiple(bool shared, TLine &bl, double by)
{
  int count = 0;
  TLine newList;

  // Point "first" at first element (may be end marker).
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last" at end marker.
  iterator last = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    iterator j = i;
    ++i; ++i;

    if (j->getSelectionFlag()) {
      if (shared && Options::warn) {
	invalidShare("MOVE");
      } else {
	// Move element and following drift to "newList".
	j->itsPosition += by;
	newList.splice(newList.end(), bl, j, i);
	++count;
	isModified = true;
      }
    } else {
      ElementBase *base = j->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("MOVE");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += moveMultiple(shared || base->isSharable(), *line, by);
      }
    }
  }

  merge(bl, newList);
  return count;
}


int Edit::moveSingleAbs(bool shared, TLine &bl, PlaceRep &it, double to)
{
  int count = 0;
  TLine newList;

  // Point "first" at first element (may be end marker).
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last" at end marker.
  iterator last = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    iterator j = i;
    ++i; ++i;
    it.enter(*j);

    if (it.isActive()) {
      if (shared && Options::warn) {
	invalidShare("MOVE");
      } else {
	// Move element and following drift to "newList".
	j->itsPosition = to;
	newList.splice(newList.end(), bl, j, i);
	++count;
	isModified = true;
      }
    } else {
      ElementBase *base = j->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("MOVE");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += moveSingleAbs(shared | base->isSharable(), *line, it, to);
      }
    }

    it.leave(*j);
  }

  merge(bl, newList);
  return count;
}


int Edit::moveSingleRel
(bool shared, TLine &bl, PlaceRep &it, PlaceRep &from, double by)
{
  int count = 0;
  TLine newList;

  // Point "first" at first element (may be end marker).
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last" at end marker.
  iterator last = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    iterator j = i;
    ++i; ++i;
    it.enter(*j);

    if (it.isActive()) {
      if (shared && Options::warn) {
	invalidShare("MOVE");
      } else {
	double fromPos = 0.0;
	from.initialize();
	for (iterator k = bl.begin(); k != bl.end(); ++k) {
	  from.enter(*k);
	  if (from.isActive()) {
	    fromPos = k->itsPosition + by;
	    goto doIt;
	  }
	  from.leave(*k);
	}
	throw OpalException("Edit::MoveSingleRel()",
			   "\"FROM\" position not found.");

      doIt:
	{
	  // Move element and following drift to "newList".
	  j->itsPosition = fromPos;
	  newList.splice(newList.end(), bl, j, i);
	  ++count;
	  isModified = true;
	}
      }
    } else {
      ElementBase *base = j->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("MOVE");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += moveSingleRel(shared || base->isSharable(),
			       *line, it, from, by);
      }
    }
    it.leave(*j);
  }

  merge(bl, newList);
  return count;
}


Edit::TLine *Edit::reflect(TLine &bl)
{
  double length = bl.back().itsPosition;
  TLine *newLine = new TLine(bl.getName());
  newLine->push_front(bl.back());

  iterator first = bl.begin();
  ++first;
  iterator last = bl.end();
  --last;
  for (iterator i = first; i != last; ++i) {
    SequenceMember member(*i);
    member.itsPosition = length - member.itsPosition;
    member.setReflectionFlag(! member.getReflectionFlag());
    newLine->push_front(member);
  }

  newLine->push_front(bl.front());
  isModified = true;
  return newLine;
}


int Edit::removeMultiple(bool shared, TLine &bl)
{
  int count = 0;

  // Point "first" at begin marker.
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last" at end marker.
  iterator last  = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    iterator j = i;
    ++i; ++i;

    if (j->getSelectionFlag()) {
      if (shared && Options::warn) {
	invalidShare("REMOVE");
	++i;
      } else {
	bl.erase(j, i);
	++count;
	isModified = true;
      }
    } else {
      ElementBase *base = j->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("REMOVE");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += removeMultiple(shared || base->isSharable(), *line);
      }
    }
  }

  return count;
}


int Edit::removeSingle(bool shared, TLine &bl, PlaceRep &it)
{
  int count = 0;

  // Point "first" at begin marker.
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last" at end marker.
  iterator last  = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    // We must get the place status before removing anything !
    it.enter(*i);
    bool doIt = it.isActive();
    it.leave(*i);

    iterator j = i;
    ++i; ++i;

    if (doIt) {
      if (shared && Options::warn) {
	invalidShare("REMOVE");
      } else {
	bl.erase(j, i);
	++count;
	isModified = true;
      }
    } else {
      ElementBase *base = j->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("REMOVE");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += removeSingle(shared || base->isSharable(), *line, it);
      }
    }
  }

  return count;
}


int Edit::replaceMultiple(bool shared, TLine &bl, ElementBase *elm)
{
  int count = 0;

  // Point "first" at first element.
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last" at end marker.
  iterator last = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    if (i->getSelectionFlag()) {
      if (shared && Options::warn) {
	invalidShare("REPLACE");
      } else {
	i->setElement(elm);

	// ada 4.5 2000 to speed up matching, add a pointer to
	// opal elements in order to avoid serching the opal elements
	i->OpalElement = Element::find(elm->getName());

	++count;
	isModified = true;
      }
    } else {
      ElementBase *base = i->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("REPLACE");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += replaceMultiple(shared || base->isSharable(), *line, elm);
      }
    }

    ++i; ++i;
  }

  return count;
}


int Edit::replaceSingle(bool shared, TLine &bl, PlaceRep &it, ElementBase *elm)
{
  int count = 0;

  // Point "first" at first element.
  iterator first = bl.begin();
  ++first; ++first;

  // Point "last" at end marker.
  iterator last = bl.end();
  --last;

  for (iterator i = first; i != last; ) {
    it.enter(*i);

    if (it.isActive()) {
      if (shared && Options::warn) {
	invalidShare("REPLACE");
      } else {
	i->setElement(elm);

	// ada 4.5 2000 to speed up matching, add a pointer to
	// opal elements in order to avoid serching the opal elements
	i->OpalElement = Element::find(elm->getName());

	++count;
	isModified = true;
      }
    } else {
      ElementBase *base = i->getElement()->removeWrappers();
      if (dynamic_cast<FlaggedBeamline *>(base)) {
	invalidLine("REPLACE");
      } else if (TLine *line = dynamic_cast<TLine *>(base)) {
	count += replaceSingle(shared || base->isSharable(), *line, it, elm);
      }
    }

    it.leave(*i);
    ++i; ++i;
  }

  return count;
}


void Edit::invalidLine(const char cmd[])
{
  if (Options::warn) {
    std::cerr << "\n### Warning ### You cannot apply \"" << cmd
	      << "\" on a nested \"LINE\".\n" << std::endl;
  }
}


void Edit::invalidShare(const char cmd[])
{
  if (Options::warn) {
    std::cerr << "\n### Warning ### You cannot apply \"" << cmd
	      << "\" on a nested shared object.\n" << std::endl;
  }
}
