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

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"

const unsigned Dim = 2;
int size = 8;
typedef Vektor<double,Dim> Vek;

Vek vabs(Vek v1)
{
  Vek v2(v1);
  for (unsigned d=0; d<Dim; ++d)
    if ( v2[d]<0 )
      v2[d] = -v2[d];
  return v2;
}

Vek vmax(Vek v1, Vek v2)
{
  Vek v0;
  for (unsigned d=0; d<Dim; ++d)
    v0[d] = max(v1[d],v2[d]);
  return v0;
}

double cutoff;
Vek vselect1(Vek v1, Vek v2)
{
  Vek v0 = vabs(v1 - v2);
  for (unsigned d=0; d<Dim; ++d)
    if ( v0[d] > cutoff )
      v0[d] = v2[d];
    else
      v0[d] = v1[d];
  return v0;
}

Vek vselect2(Vek diff, Vek cutoff)
{
  Vek v0;
  for (unsigned d=0; d<Dim; ++d)
    if ( cutoff[d] >= fabs(diff[d])  )
      v0[d] = diff[d];
    else
      v0[d] = 0.0;
  return v0;
}

double el0(Vek v1)
{
  return v1[0];
}

UNARY_FUNCTION(double,el0,Vek)
UNARY_FUNCTION(Vek,vabs,Vek)
BINARY_FUNCTION(Vek,vmax,Vek, Vek)
BINARY_FUNCTION(Vek,vselect1,Vek,Vek)
BINARY_FUNCTION(Vek,vselect2,Vek,Vek)

int main(int argc, char *argv[])
{
  // initialize Ippl, and create Inform object for output messages
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  Index I(size), J(size);
  int i,j;
  FieldLayout<Dim> layout(I,J);
#ifdef __MWERKS__
  typedef Field< Vek ,Dim,UniformCartesian<Dim>, 
    UniformCartesian<Dim>::DefaultCentering> FV;
  typedef Field< double ,Dim,UniformCartesian<Dim>,
    UniformCartesian<Dim>::DefaultCentering> Fd;
#else
  typedef Field< Vek ,Dim,UniformCartesian<Dim> > FV;
  typedef Field< double ,Dim,UniformCartesian<Dim> > Fd;
#endif // __MWERKS__
  FV V1(layout),V2(layout),V3(layout);
  Fd d1(layout),d2(layout),d3(layout);
  double err;
  
  for (i=0; i<size; ++i)
    for (j=0; j<size; ++j) {
      V1[i][j] = Vek( double(i), double(j) ) ;
      V2[i][j] = Vek( double(j), double(i) ) ;
    }

  V3 = 0.0;
  assign( V3, vabs(V2-V1)+V3+vabs(V3) );

  err = 0;
  for (i=0; i<size; ++i)
    for (j=0; j<size; ++j)
      {
	Vek v = V3[i][j].get();
	v -= Vek( abs(j-i), abs(i-j) );
	err += fabs(v[0]);
	err += fabs(v[1]);
      }
  if ( err==0 )
    testmsg << "PASSED UNARY" << endl;
  else
    testmsg << "FAILED UNARY" << endl;
  
  assign( V3, vmax(V2,V1) );
  err = 0;
  for (i=0; i<size; ++i)
    for (j=0; j<size; ++j)
      {
	Vek v = V3[i][j].get();
	double m = (i>j) ? i : j;
	v -= Vek( m,m );
	err += fabs(v[0]);
	err += fabs(v[1]);
      }
  if ( err==0 )
    testmsg << "PASSED BINARY" << endl;
  else
    testmsg << "FAILED BINARY" << endl;

  cutoff = 1.0;
  assign( V3, vselect1(V1,V2) );
  err = 0;
  for (i=0; i<size; ++i)
    for (j=0; j<size; ++j)
      {
	Vek v = V3[i][j].get();
	v -= vselect1( Vek(i,j), Vek(j,i) );
	err += fabs(v[0]);
	err += fabs(v[1]);
      }
  if ( err==0 )
    testmsg << "PASSED SELECT1" << endl;
  else
    testmsg << "FAILED SELECT1" << endl;

  V3 = 1.0;
  V2 += vselect2(V1-V2,V3) ;
  err = 0;
  for (i=0; i<size; ++i)
    for (j=0; j<size; ++j)
      {
	Vek v = V2[i][j].get();
	v -= vselect1( Vek(i,j), Vek(j,i) );
	err += fabs(v[0]);
	err += fabs(v[1]);
      }
  if ( err==0 )
    testmsg << "PASSED SELECT2" << endl;
  else
    testmsg << "FAILED SELECT2" << endl;

  assign( d1, fabs(el0(V2)-el0(V1)) );

  return 0;
}
//----------------------------------------------------------------------
/***************************************************************************
 * $RCSfile: tz.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:42 $
 * IPPL_VERSION_ID: $Id: tz.cpp,v 1.1.1.1 2000/11/30 21:06:42 adelmann Exp $ 
 ***************************************************************************/
