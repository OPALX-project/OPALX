// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Utility/PAssert.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/BareField.h"
#include "Field/GuardCellSizes.h"


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg;
  // testmsg.setPrintNode();

  const int Dim=2;
  Index I(10),J(10);
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL,4);
  GuardCellSizes<Dim> gc(1);
  typedef BareField<int,Dim> F;
  F A(layout,gc);

  // A should be constructed compressed.
  bool passed_this;
  F::iterator_if lf;
  int count;

  //////////////////////////////////////////////////////////////////////
  // Test if it is constructed compressed 
  // (or uncompressed, if --nofieldcompression)
  //////////////////////////////////////////////////////////////////////
  passed_this = true;

  if (IpplInfo::noFieldCompression) {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) 
      {
	passed_this = passed_this && !(*lf).second->IsCompressed();
	if ( (*lf).second->IsCompressed() )
	  testmsg << "FAILED: An LField is compressed," << count << endl;
      }
    if ( passed_this )
      testmsg << "PASSED: Field is constructed uncompressed, " << count 
	      << endl;
  }
  else {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) 
      {
	passed_this = passed_this && (*lf).second->IsCompressed();
	if ( ! (*lf).second->IsCompressed() )
	  testmsg << "FAILED: An LField is uncompressed," << count << endl;
      }
    if ( passed_this )
      testmsg << "PASSED: Field is constructed compressed, " << count << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether fillGuardCells destroys the compression.
  //////////////////////////////////////////////////////////////////////
  passed_this = true;
  A.fillGuardCells();
  if (IpplInfo::noFieldCompression) {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count)
      {
	passed_this = passed_this && !(*lf).second->IsCompressed();
	if ( (*lf).second->IsCompressed() )
	  testmsg << "FAILED: Compressed after fillGuardCells, " << count 
		  << endl;
      }
    if ( passed_this )
      testmsg << "PASSED: Uncompressed after fillGuardCells, " << count
	      << endl;
  }
  else {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count)
      {
	passed_this = passed_this && (*lf).second->IsCompressed();
	if ( ! (*lf).second->IsCompressed() )
	  testmsg << "FAILED: Uncompressed after fillGuardCells, " << count 
		  << endl;
      }
    if ( passed_this )
      testmsg << "PASSED: Compressed after fillGuardCells, " << count << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether assigning an index uncompresses.
  //////////////////////////////////////////////////////////////////////
  passed_this = true;
  assign(A[I][J] , I + 10*J);
  for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count)
    {
      passed_this = passed_this && !(*lf).second->IsCompressed();
      if ( (*lf).second->IsCompressed() )
	testmsg << "FAILED: Compressed after assigning Index, " << count 
		<< endl;
    }
  if ( passed_this )
    testmsg << "PASSED: Uncompressed after assigning Index, " << count
	    << endl;
  if (IpplInfo::noFieldCompression) {
    {
      int s = sum(A);
      int il = I.length();
      int jl = J.length();
      int ss = jl*il*(il-1)/2 + il*jl*(jl-1)*5;
      if ( s != ss )
	{
	  testmsg << "FAILED: incorrect sum." << endl;
	  testmsg << "A=" << A << endl;
	  testmsg << "sum(A) = " << s << endl;
	  testmsg << "calc = " << ss << endl;
	}
      else 
	testmsg << "PASSED: correct sum." << endl;
    }
  }
  else {
    {
      int s = sum(A);
      int il = I.length();
      int jl = J.length();
      int ss = jl*il*(il-1)/2 + il*jl*(jl-1)*5;
      if ( s != ss )
	{
	  testmsg << "FAILED: incorrect sum." << endl;
	  testmsg << "A=" << A << endl;
	  testmsg << "sum(A) = " << s << endl;
	  testmsg << "calc = " << ss << endl;
	}
      else 
	testmsg << "PASSED: correct sum." << endl;
    }
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether assigning a constant compresses.
  //////////////////////////////////////////////////////////////////////
  passed_this = true;
  A = 1;
  if (IpplInfo::noFieldCompression) {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
      passed_this = passed_this && !(*lf).second->IsCompressed();
      if ( (*lf).second->IsCompressed() )
	testmsg << "FAILED: Compressed after assigning constant, " << count 
		<< endl;
    }
    if ( passed_this )
      testmsg << "PASSED: Uncompressed after assigning constant, " << count 
	      << endl;
  }
  else {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count) {
      passed_this = passed_this && (*lf).second->IsCompressed();
      if ( !(*lf).second->IsCompressed() )
	testmsg << "FAILED: Uncompressed after assigning constant, " << count 
		<< endl;
    }
    if ( passed_this )
      testmsg << "PASSED: Compressed after assigning constant, " << count 
	      << endl;
  }
  {
    int s = sum(A);
    if ( s != I.length()*J.length() )
      testmsg << "FAILED: incorrect sum." << endl;
    else 
      testmsg << "PASSED: correct sum." << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether assigning a constant with indexes compresses.
  //////////////////////////////////////////////////////////////////////
  passed_this = true;
  // First uncompress it.
  assign(A[I][J] , I+10*J);
  // Then compress it.
  assign(A[I][J] , 1 );
  if (IpplInfo::noFieldCompression) {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count)
      {
	passed_this = passed_this && !(*lf).second->IsCompressed();
	if ( (*lf).second->IsCompressed() )
	  testmsg << "FAILED: Compressed after I-assigning constant, " 
		  << count << endl;
      }
    if ( passed_this )
      testmsg << "PASSED: Unompressed after I-assigning constant, " << count 
	      << endl;
  }
  else {
    for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count)
      {
	passed_this = passed_this && (*lf).second->IsCompressed();
	if ( !(*lf).second->IsCompressed() )
	  testmsg << "FAILED: Uncompressed after I-assigning constant, " 
		  << count << endl;
      }
    if ( passed_this )
      testmsg << "PASSED: Compressed after I-assigning constant, " << count 
	      << endl;
  }
  {
    int s = sum(A);
    if ( s != I.length()*J.length() )
      testmsg << "FAILED: incorrect sum." << endl;
    else 
      testmsg << "PASSED: correct sum." << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether assigning a subrange uncompresses.
  //////////////////////////////////////////////////////////////////////
  passed_this = true;
  Index I1( I.min()+1 , I.max()-1 );
  Index J1( J.min()+1 , J.max()-1 );
  A=0;
  assign(A[I1][J1] , 1 );
  for (lf=A.begin_if(), count=0; lf!=A.end_if(); ++lf, ++count)
    {
      passed_this = passed_this && !(*lf).second->IsCompressed();
      if ( (*lf).second->IsCompressed() )
	testmsg << "FAILED: Compressed after assigning subrange, " << count 
		<< endl;
    }
  if ( passed_this )
    testmsg << "PASSED: Uncompressed after assigning subrange, " << count 
	    << endl;
  {
    int s = sum(A);
    if ( s != I1.length()*J1.length() )
      testmsg << "FAILED: incorrect sum." << endl;
    else 
      testmsg << "PASSED: correct sum." << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether assigning a single element uncompresses one.
  //////////////////////////////////////////////////////////////////////
  A = 0;
  if (!IpplInfo::noFieldCompression) {
    for (lf=A.begin_if(); lf!=A.end_if(); ++lf)
      PAssert( (*lf).second->IsCompressed() );
  }
  assign(A[3][3] , 1 );
  count = 0;
  int reduced_count = 0;
  if (IpplInfo::noFieldCompression) {
    if (A.CompressedFraction() != 0) {
      testmsg << "FAILED: Compressed somewhere after assigning single"
	      << endl;
    }
    else {
      testmsg << "PASSED: Still uncompressed after assigning single" << endl;
    }
  }
  else {
    for (lf=A.begin_if(); lf!=A.end_if(); ++lf)
      {
	if ( (*lf).second->IsCompressed() )
	  ++count;
      }
    reduced_count = 0;
    reduce(&count,&count+1,&reduced_count,OpAddAssign());
    if ( reduced_count == 3 )
      testmsg << "PASSED: Uncompressed one after assigning single" << endl;
    else
      testmsg << "FAILED: " << count << " compressed after assigning single." 
	      << endl;
  }
  {
    int s = sum(A);
    if ( s != 1 )
      testmsg << "FAILED: incorrect sum." << endl;
    else 
      testmsg << "PASSED: correct sum." << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether assigning a single element uncompresses one.
  //////////////////////////////////////////////////////////////////////
  A = 0;
  if (!IpplInfo::noFieldCompression) {
    for (lf=A.begin_if(); lf!=A.end_if(); ++lf)
      PAssert( (*lf).second->IsCompressed() );
  }
  assign(A[4][3] , 1 );
  count = 0;
  if (IpplInfo::noFieldCompression) {
    if (A.CompressedFraction() != 0) {
      testmsg << "FAILED: Compressed somewhere after assigning single"
	      << endl;
    }
    else {
      testmsg << "PASSED: Still uncompressed after assigning single" << endl;
    }
  }
  else {
    for (lf=A.begin_if(); lf!=A.end_if(); ++lf)
      {
	if ( (*lf).second->IsCompressed() )
	  ++count;
      }
    reduced_count = 0;
    reduce(&count,&count+1,&reduced_count,OpAddAssign());
    if (Ippl::deferGuardCellFills) {
      if ( reduced_count == 3 )
	testmsg << "PASSED: Uncompressed one after assigning in guard cell"
		<< endl;
      else
	testmsg << "FAILED: " << count
		<< " compressed after assigning in guard cell." << endl;
    }
    else {
      if ( reduced_count == 2 )
	testmsg << "PASSED: Uncompressed two after assigning in guard cell"
		<< endl;
      else
	testmsg << "FAILED: " << count
		<< " compressed after assigning in guard cell." << endl;
    }
  }
  {
    int s = sum(A);
    if ( s != 1 )
      testmsg << "FAILED: incorrect sum." << endl;
    else 
      testmsg << "PASSED: correct sum." << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;

  //////////////////////////////////////////////////////////////////////
  // Test whether an operation on that array leaves it correct.
  //////////////////////////////////////////////////////////////////////
  A *= A;
  count = 0;
  if (IpplInfo::noFieldCompression) {
    if (A.CompressedFraction() != 0) {
      testmsg << "FAILED: Compressed somewhere after squaring" << endl;
    }
    else {
      testmsg << "PASSED: Still uncompressed after squaring" << endl;
    }
  }
  else {
    for (lf=A.begin_if(); lf!=A.end_if(); ++lf)
      {
	if ( (*lf).second->IsCompressed() )
	  ++count;
      }
    reduced_count = 0;
    reduce(&count,&count+1,&reduced_count,OpAddAssign());
    if (Ippl::deferGuardCellFills) {
      if ( reduced_count == 3 )
	testmsg << "PASSED: Uncompressed one after squaring"
		<< endl;
      else
	testmsg << "FAILED: " << count
		<< " compressed after squaring." << endl;
    }
    else {
      if ( reduced_count == 2 )
	testmsg << "PASSED: Uncompressed two after squaring"
		<< endl;
      else
	testmsg << "FAILED: " << count
		<< " compressed after squaring." << endl;
    }
  }
  {
    int s = sum(A);
    if ( s != 1 )
      testmsg << "FAILED: incorrect sum." << endl;
    else 
      testmsg << "PASSED: correct sum." << endl;
  }
  testmsg << "        compressed fraction=" << A.CompressedFraction() << endl;
  
  //////////////////////////////////////////////////////////////////////
  // Make sure we can construct a Field of Maps.  
  //////////////////////////////////////////////////////////////////////
  BareField< map<int,double> , Dim > B(layout);
  BareField< map<int,double> , Dim >::iterator_if lb;
  map<int,double> m;
  m[1] = cos(1.0);
  m[2] = cos(2.0);
  B[2][2] = m;
  count = 0;
  if (IpplInfo::noFieldCompression) {
    if (B.CompressedFraction() != 0) {
      testmsg << "FAILED: Compressed somewhere in Field of maps" << endl;
    }
    else {
      testmsg << "PASSED: Still uncompressed Field of maps" << endl;
    }
  }
  else {
    for (lb=B.begin_if(); lb!=B.end_if(); ++lb)
      {
	if ( (*lb).second->IsCompressed() )
	  ++count;
      }
    reduced_count = 0;
    reduce(&count,&count+1,&reduced_count,OpAddAssign());
    if ( reduced_count == 3 )
      testmsg << "PASSED: Uncompressed one in Field of maps" << endl;
    else
      testmsg << "FAILED: " <<  count
	      << " compressed after in Field of maps." << endl;
  }
  testmsg << "        compressed fraction=" << B.CompressedFraction() << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: compressed.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: compressed.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/

