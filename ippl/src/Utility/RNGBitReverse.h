// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef RNG_BIT_REVERSE_H
#define RNG_BIT_REVERSE_H

/***********************************************************************
 * 
 * class RNGBitReverse
 * class RNGBitReverseSequence : public SequenceGen<RNGBitReverse>
 *
 * RNGBitReverse generates a set of values distributed over some domain
 * based on a bit-reversal algorithm.  These numbers will be uniformly
 * distributed over the [0...1] domain.  It requires the following
 * parameters and values:
 *
 * base - used in bit reversal calculation; must be a prime number
 * seed - a starting unsigned integer, which is incremented as new values
 *        are calculated.  If not specified, it defaults to 1 initially.
 *        It is incremented after each calculation.
 *
 * RNGBitReverseSequence is derived from SequenceGen, and makes it easier
 * to use this RNG in expressions.
 * Use RNGBitReverse as a scalar or container element, and use
 * RNGBitReverseSequence when you need a sequence of numbers to fill
 * a container.
 *
 ***********************************************************************/

// include files
#include "Utility/SequenceGen.h"


class RNGBitReverse {

public:
  // return type
  typedef double Return_t;

public:
  // default constructor
  RNGBitReverse(unsigned int base = 2, unsigned long seed = 1)
    : Base(base), Seed(seed) {}

  // copy constructor
  RNGBitReverse(const RNGBitReverse& brg)
    : Base(brg.Base), Seed(brg.Seed) {}

  // destructor
  ~RNGBitReverse(void) {}

  // advance seed for RNG n times
  inline void AdvanceSeed(unsigned long n = 0) { Seed += n; }

  // set seed to specified value
  inline void SetSeed(unsigned long seed)        { Seed = seed; }

  // set base to specified value (should be a prime number!)
  inline void SetBase(unsigned int base)         { Base = base; }

  // get seed value
  inline unsigned long GetSeed(void) const       { return Seed; }

  // get base value
  inline unsigned int GetBase(void) const        { return Base; }

  // get the next random number in sequence
  inline Return_t GetRandom(void) const {
    Return_t rev = 0.0;
    Return_t power = 1.0;
    long inum = Seed;
    while (inum > 0) {
      int iquot = inum/Base;
      int irem = inum - iquot*Base;
      power /= Return_t(Base);
      rev += Return_t(irem)*power;
      inum = iquot;
    }
    Seed++;
    return rev;
  }

  // pseudonym for GetRandom()
  inline Return_t operator()(void) const { return GetRandom(); }

  // conversion to Return_t, same as GetRandom()
  inline operator Return_t() const { return GetRandom(); }

private:
  // the base and seed for this generator
  unsigned int          Base;
  mutable unsigned long Seed;
};

RNG_BASIC_MATH(RNGBitReverse)


// A version of SequenceGen with extra constructors to make using this
// class easier.  This is the version that people should use to fill
// containers with a random number sequence in an expression.  This
// class is PETE-aware via its inheritance from SequenceGen.

class RNGBitReverseSequence : public SequenceGen<RNGBitReverse> {

public:
  // default constructor
  RNGBitReverseSequence(unsigned int base = 2, unsigned long seed = 1)
    : SequenceGen<RNGBitReverse>(RNGBitReverse(base,seed)) {}

  // copy constructor
  RNGBitReverseSequence(const RNGBitReverseSequence& rngseq)
    : SequenceGen<RNGBitReverse>(rngseq.getGenerator()) {}

  // destructor
  ~RNGBitReverseSequence(void) {}

  // wrappers around RNG generator functions
  inline void AdvanceSeed(unsigned long adv = 0) {
    getGenerator().AdvanceSeed(adv);
  }

  inline void SetSeed(unsigned long seed) { getGenerator().SetSeed(seed); }
  inline void SetBase(unsigned int base)  { getGenerator().SetBase(base); }
  inline unsigned long GetSeed(void) const { return getGenerator().GetSeed(); }
  inline unsigned int GetBase(void) const { return getGenerator().GetBase(); }

  inline Return_t GetRandom(void) { return getGenerator().GetRandom(); }
  inline Return_t operator()(void) { return getGenerator().GetRandom(); }
};

#endif // RNG_BIT_REVERSE_H

/***************************************************************************
 * $RCSfile: RNGBitReverse.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: RNGBitReverse.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
