// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef RNG_XCI_H
#define RNG_XCI_H

/***********************************************************************
 * 
 * class RNGXCI
 * class RNGXCISequence : public SequenceGen<RNGXCI>
 *
 * RNGXCI is a simple class that implements random number generator from LANL
 * XCI group (Forrest Brown), in the range (0...1).  The first class is used 
 * as a scalar, as an element of a Field or ParticleAttrib, or as a template
 * parameter to SequenceGen.  The second class is derived from SequenceGen,
 * and makes it easier to use this RNG in an expression as a sequence of
 * numbers that will fill a Field or ParticleAttrib container.
 * Use RNGXCI as a scalar or container element, and use
 * RNGXCISequence when you need a sequence of numbers to fill a container.
 *
 ***********************************************************************/

// include files
#include "Utility/SequenceGen.h"

#ifdef IPPL_USE_SINGLE_PRECISION
 #define CTYPE float
#else
 #define CTYPE double
#endif

// define type that we ensure is 8 bytes long
#define LONG_IS_8_BYTES (((1L<<16)<<16)<<16)
typedef long long  RNlong;


class RNGXCI {

public:
  // return type
  typedef CTYPE Return_t;

public:
  // default constructor
  RNGXCI(RNlong advance = 0L)
    : Seed(1L) { AdvanceSeed(advance); }

  // copy constructor
  RNGXCI(const RNGXCI& rng)
    : Seed(rng.Seed) {}

  // destructor
  ~RNGXCI(void) {}

  // advance seed for random number generator n times
  inline void AdvanceSeed(RNlong n = 0L) {
    //  skip ahead n RNs:   Seed*RN_MULT^n mod RN_MOD
    while( n<0 )  n += RN_PERIOD;  // add period till >0
    n &= RN_MASK;                  // mod RN_MOD
    RNlong  gen=1L,  g=RN_MULT;    // get gen=RN_MULT^n,  in
    for( ; n; n>>=1 ) {            //    log2(n) ops, not n ops !
      if( n&1 ) gen = gen*g & RN_MASK;
      g = g*g & RN_MASK;
    }
    Seed = gen*Seed & RN_MASK;
    return;
  }

  // set seed to user-specified value
  inline void SetSeed(RNlong s) {
    Seed = (s|1) & RN_MASK;        // must be odd to get full period
  }

  // get seed value
  inline RNlong GetSeed(void) const { return Seed; }

  // return the next pseudo-random number (from 0 ... 1)
  inline Return_t GetRandom(void) const {
    Return_t r = Seed*RN_NORM;
    Seed = Seed*RN_MULT & RN_MASK;
    return r;
  }

  // pseudonym for GetRandom()
  inline Return_t operator()(void) const { return GetRandom(); }

  // conversion to Return_t, same as GetRandom()
  inline operator Return_t() const { return GetRandom(); }

  // return the period of the RNG
  static Return_t GetRandMax(void) { return Return_t(RN_MASK); }

private:
  mutable RNlong Seed;

  // constants for RN gen:       Seed = Seed*RN_MULT mod RN_MOD
  static const RNlong RN_MULT;   // 5^19 
  static const RNlong RN_MOD;    // 2^48 
  static const RNlong RN_PERIOD; // period
  static const double RN_NORM;   // normalize to (0,1)
  static const RNlong RN_MASK;   // 48-bit mask
};

RNG_BASIC_MATH(RNGXCI)


// A version of SequenceGen with extra constructors to make using this
// class easier.  This is the version that people should use to fill
// containers with a random number sequence in an expression.  This
// class is PETE-aware via its inheritance from SequenceGen.

class RNGXCISequence : public SequenceGen<RNGXCI> {

public:
  // default constructor
  RNGXCISequence(int advance = 0)
    : SequenceGen<RNGXCI>(RNGXCI(advance)) {}

  // copy constructor
  RNGXCISequence(const RNGXCISequence& rngseq)
    : SequenceGen<RNGXCI>(rngseq.getGenerator()) {}

  // destructor
  ~RNGXCISequence(void) {}

  // wrappers around RNG generator functions
  inline void AdvanceSeed(RNlong adv = 0L) {
    getGenerator().AdvanceSeed(adv);
  }
  inline void     SetSeed(RNlong seed) { getGenerator().SetSeed(seed); }
  inline RNlong GetSeed(void) const { return getGenerator().GetSeed(); }
  inline Return_t GetRandom(void) { return getGenerator().GetRandom(); }
  inline Return_t operator()(void) { return getGenerator().GetRandom(); }
  static Return_t GetRandMax(void) { return RNGXCI::GetRandMax(); }
};


#endif // RNG_XCI_H

/***************************************************************************
 * $RCSfile: RNGXCI.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: RNGXCI.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/

