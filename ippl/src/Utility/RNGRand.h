// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef RNG_RAND_H
#define RNG_RAND_H

/***********************************************************************
 * 
 * class RNGRand
 * class RNGRandSequence : public SequenceGen<RNGRand>
 *
 * Simple wrapper around C rand() function, to generate random numbers
 * in the range [0...1].  The first class may be
 * used standalone, or as a template parameter to SequenceGen.  The
 * second class is derived from SequenceGen, and makes it easier to
 * use this RNG in expressions.
 * Use RNGRand as a scalar or container element, and use
 * RNGRandSequence when you need a sequence of numbers to fill a container.
 *
 ***********************************************************************/

// include files
#include "Utility/SequenceGen.h"

#ifdef IPPL_USE_SINGLE_PRECISION
 #define CTYPE float
#else
 #define CTYPE double
#endif

#include <cstdlib>

class RNGRand {

public:
  // return type
  typedef CTYPE Return_t;

public:
  // default constructor
  RNGRand(int advance = 0) : CurrentSeed(1U) { AdvanceSeed(advance); }

  // copy constructor
  RNGRand(const RNGRand& rng) : CurrentSeed(rng.CurrentSeed) {}

  // destructor
  ~RNGRand(void) {}

  // advance indicates number of times to advance random number source
  // there's no standard way to do this with rand(), so make one up!
  // just require one-to-one correspondence of old and new seeds
  inline void AdvanceSeed(int advance = 0) {
    for (int iadv=0; iadv<advance; iadv++) {
      // first reseed with current seed
      srand(CurrentSeed);
      // use first random number as seed step
      int seedstep = rand();
      CurrentSeed = (CurrentSeed + seedstep) % RAND_MAX;
    }
    // finally, reseed with new seed value
    srand(CurrentSeed);
  }

  // set seed to user-specified value (any shifting is done by srand)
  inline void SetSeed(unsigned int seed) {
    CurrentSeed = seed;
    srand(CurrentSeed);
  }

  // get seed value
  inline unsigned int GetSeed(void) const { return CurrentSeed; }

  // return the next pseudo-random number (from 0 ... 1)
  inline Return_t GetRandom(void) const {
    return ( (Return_t(rand())) / (Return_t(RAND_MAX)+1) );
  }

  // pseudonym for GetRandom()
  inline Return_t operator()(void) const { return GetRandom(); }

  // conversion to Return_t, same as GetRandom()
  inline operator Return_t() const { return GetRandom(); }

  // return the period of the RNG
  static Return_t GetRandMax(void) { return Return_t(RAND_MAX)+1; }

private:
  unsigned int CurrentSeed;
};

RNG_BASIC_MATH(RNGRand)


// A version of SequenceGen with extra constructors to make using this
// class easier.  This is the version that people should use to fill
// containers with a random number sequence in an expression.  This
// class is PETE-aware via its inheritance from SequenceGen.

class RNGRandSequence : public SequenceGen<RNGRand> {

public:
  // default constructor
  RNGRandSequence(int advance = 0)
    : SequenceGen<RNGRand>(RNGRand(advance)) {}

  // copy constructor
  RNGRandSequence(const RNGRandSequence& rngseq)
    : SequenceGen<RNGRand>(rngseq.getGenerator()) {}

  // destructor
  ~RNGRandSequence(void) {}

  // wrappers around RNG generator functions
  inline void     AdvanceSeed(int adv = 0) { getGenerator().AdvanceSeed(adv); }
  inline void     SetSeed(unsigned int seed) { getGenerator().SetSeed(seed); }
  inline unsigned int GetSeed(void) const { return getGenerator().GetSeed(); }
  inline Return_t GetRandom(void) { return getGenerator().GetRandom(); }
  inline Return_t operator()(void) { return getGenerator().GetRandom(); }
  static Return_t GetRandMax(void) { return RNGRand::GetRandMax(); }
};


#endif // RNG_RAND_H

/***************************************************************************
 * $RCSfile: RNGRand.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: RNGRand.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/

