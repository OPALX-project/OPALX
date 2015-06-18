// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef PASSERT_H
#define PASSERT_H

#include <exception>
#include <stdexcept>
//////////////////////////////////////////////////////////////////////
//
// This is a compile time assert.
// That is, if you say:
//   CTAssert<true>::test();
// it compiles just fine and inserts no code.
// If you say:
//   CTAssert<false>::test();
// you get a compile error that it can't find CTAssert<false>::test().
//
// The template argument can of course be a calculation of const bools
// that are known at compile time.
//
//////////////////////////////////////////////////////////////////////

template<bool B> struct IpplCTAssert {};

template<> struct IpplCTAssert<true> { static void test() {} };

#if defined(NOCTAssert)
#define CTAssert(c)
#else
#define CTAssert(c) IpplCTAssert<(c)>::test()
#endif

//===========================================================================//
// class assertion - exception notification class for assertions

// This class should really be derived from std::runtime_error, but
// unfortunately we don't have good implementation of the library standard
// yet, on compilers other than KCC.  So, this class will keep with the
// "what" method evidenced in the standard, but dispense with inheriting from
// classes for which we don't have implementations...
//===========================================================================//

class assertion: public std::runtime_error
{
    char *msg;
public:
    assertion( const char *cond, const char *file, int line );
    assertion( const char *m );
    assertion( const assertion& a );
  ~assertion() throw() { delete[] msg; }
    assertion& operator=( const assertion& a );
    const char* what() const noexcept { return msg; }
};

//---------------------------------------------------------------------------//
// Now we define a run time assertion mechanism.  We will call it "PAssert",
// to reflect the idea that this is for use in IPPL per se, recognizing that
// there are numerous other assertion facilities in use in client codes.
//---------------------------------------------------------------------------//

// These are the functions that will be called in the assert macros.

void toss_cookies( const char *cond, const char *file, int line );
void insist( const char *cond, const char *msg, const char *file, int line );

//---------------------------------------------------------------------------//
// The PAssert macro is intended to be used for validating preconditions
// which must be true in order for following code to be correct, etc.  For
// example, PAssert( x > 0. ); y = sqrt(x);  If the assertion fails, the code
// should just bomb.  Philosophically, it should be used to feret out bugs in
// preceding code, making sure that prior results are within reasonable
// bounds before proceeding to use those results in further computation, etc.
//---------------------------------------------------------------------------//

#ifdef NOPAssert
#define PAssert(c)
#else
#define PAssert(c) if (!(c)) toss_cookies( #c, __FILE__, __LINE__ );
#endif

//---------------------------------------------------------------------------//
// The PInsist macro is akin to the PAssert macro, but it provides the
// opportunity to specify an instructive message.  The idea here is that you
// should use Insist for checking things which are more or less under user
// control.  If the user makes a poor choice, we "insist" that it be
// corrected, providing a corrective hint.
//---------------------------------------------------------------------------//

#define PInsist(c,m) if (!(c)) insist( #c, m, __FILE__, __LINE__ );

//---------------------------------------------------------------------------//
// NOTE:  We provide a way to eliminate assertions, but not insistings.  The
// idea is that PAssert is used to perform sanity checks during program
// development, which you might want to eliminate during production runs for
// performance sake.  PInsist is used for things which really really must be
// true, such as "the file must've been opened", etc.  So, use PAssert for
// things which you want taken out of production codes (like, the check might
// inhibit inlining or something like that), but use PInsist for those things
// you want checked even in a production code.
//---------------------------------------------------------------------------//

#endif // PASSERT_H

/***************************************************************************
 * $RCSfile: PAssert.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: PAssert.h,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $
 ***************************************************************************/
