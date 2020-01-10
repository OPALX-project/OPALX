// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 ***************************************************************************/

#ifndef RNG_SIMPLE_H
#define RNG_SIMPLE_H

/***********************************************************************
 * 
 * class RNGSimple
 * class RNGSimpleSequence : public SequenceGen<RNGSimple>
 *
 * Simple class that implements random number generator a la
 * Numerical Recipes, in the range [0...1].  The first class may be
 * used standalone, or as a template parameter to SequenceGen.  The
 * second class is derived from SequenceGen, and makes it easier to
 * use this RNG in expressions.
 * Use RNGSimple as a scalar or container element, and use
 * RNGSimpleSequence when you need a sequence of numbers to fill a container.
 *
 ***********************************************************************/

#include "Utility/SequenceGen.h"

class RNGSimple {

public:
  typedef double Return_t;

public:
  // default constructor
  RNGSimple(int advance = 0)
    : CurrentSeed(RandShift + 1), CurrentRand(0) { AdvanceSeed(advance); }

  // copy constructor
  RNGSimple(const RNGSimple& rng)
    : CurrentSeed(rng.CurrentSeed), CurrentRand(rng.CurrentRand) { }

  // destructor
  ~RNGSimple(void) {}

  //   advance indicates number of times to advance random number source
  inline void AdvanceSeed(int advance = 0) {
    for (int iadv=0; iadv<advance; iadv++) 
      CurrentSeed = (CurrentSeed + RandShift) % RandModulus;
    CurrentRand = CurrentSeed;
  }

  // set seed to user-specified value, plus shift to ensure it is large
  inline void SetSeed(unsigned long seed) {
    CurrentSeed = (seed + RandShift) % RandModulus;
    CurrentRand = CurrentSeed;
  }

  // get seed value
  inline unsigned long GetSeed(void) const { return CurrentSeed; }

  // return the next pseudo-random number (from 0 ... 1)
  inline Return_t GetRandom(void) const {
    CurrentRand = (CurrentRand * RandMultipplier + RandShift) % RandModulus;
    return ( Return_t(CurrentRand) / Return_t(RandModulus) );
  }

  // pseudonym for GetRandom()
  inline Return_t operator()(void) const { return GetRandom(); }

  // conversion to Return_t, same as GetRandom()
  inline operator Return_t() const { return GetRandom(); }

  // return the period of the RNG
  static Return_t GetRandMax(void) { return Return_t(RandModulus); }

private:
  long CurrentSeed;
  mutable long CurrentRand;

  static const long RandModulus;
  static const long RandMultipplier;
  static const long RandShift;
};

RNG_BASIC_MATH(RNGSimple)


// A version of SequenceGen with extra constructors to make using this
// class easier.  This is the version that people should use to fill
// containers with a random number sequence in an expression.  This
// class is PETE-aware via its inheritance from SequenceGen.

class RNGSimpleSequence : public SequenceGen<RNGSimple> {

public:
  // default constructor
  RNGSimpleSequence(int advance = 0)
    : SequenceGen<RNGSimple>(RNGSimple(advance)) {}

  // copy constructor
  RNGSimpleSequence(const RNGSimpleSequence& rngseq)
    : SequenceGen<RNGSimple>(rngseq.getGenerator()) {}

  // destructor
  ~RNGSimpleSequence(void) {}

  // wrappers around RNG generator functions
  inline void     AdvanceSeed(int adv = 0) { getGenerator().AdvanceSeed(adv); }
  inline void     SetSeed(unsigned long seed) { getGenerator().SetSeed(seed); }
  inline unsigned long GetSeed(void) const { return getGenerator().GetSeed(); }
  inline Return_t GetRandom(void) { return getGenerator().GetRandom(); }
  inline Return_t operator()(void) { return getGenerator().GetRandom(); }
  static Return_t GetRandMax(void) { return RNGSimple::GetRandMax(); }
};


#endif // RNG_SIMPLE_H

