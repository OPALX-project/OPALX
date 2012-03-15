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

// TestFieldDebug2.cpp , Tim Williams 10/08/1997
// This tests the revised functions [gg]fp3(), which internally tests the
// revised sfp3(), which are all supposed to allow printing including
// global guard layers, even with multipple vnodes.  Also tests the
// vectorFace-centered-appropriate BC application in 3D, with multipple vnodes.
// Not a self-checking test; requires tedious exmaination of output, but it's
// worth it if there is a suspected problem with any of this because of changes
// in IPPL.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "Field/BCond.h"
#include "Field/GuardCellSizes.h"
#include "AppTypes/Vektor.h"
#include "Utility/FieldDebug.h"

//-----------------------------------------------------------------------------
// User-inserted prototypes to get debugger access (examples for user ref.):
// Scalar (double):
void dfp3(BareField<double,3U> f) {fp3(f);}
void defp3(BareField<double,3U> f, int i, int j, int k) {efp3(f,i,j,k);}
void dsfp3(BareField<double,3U>& f,
	   int base1, int bound1, int stride1,
	   int base2, int bound2, int stride2,
	   int base3, int bound3, int stride3) {
  sfp3(f,base1,bound1,stride1,base2,bound2,stride2,base3,bound3,stride3);}
// Vektor (double):
void vdfp3(BareField<Vektor<double,3U>,3U> f) {fp3(f);}
void vdefp3(BareField<Vektor<double,3U>,3U> f, int i, int j, int k) {
  efp3(f,i,j,k);}
void vdsfp3(BareField<Vektor<double,3U>,3U>& f,
	    int base1, int bound1, int stride1,
	    int base2, int bound2, int stride2,
	    int base3, int bound3, int stride3) {
  sfp3(f,base1,bound1,stride1,base2,bound2,stride2,base3,bound3,stride3);}
//-----------------------------------------------------------------------------


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  const unsigned D=3U;
  int nvnodes = 4;
  int d, face;
  unsigned nverts[D]; for (d=0; d<D; d++) nverts[d] = 5;
  unsigned ncells[D]; for (d=0; d<D; d++) ncells[d] = nverts[d] - 1;
  unsigned nguards[D]; for (d=0; d<D; d++) nguards[d] = 2;

  NDIndex<D> verts; for (d=0; d<D; d++) verts[d] = Index(nverts[d]);
  NDIndex<D> cells; for (d=0; d<D; d++) cells[d] = Index(ncells[d]);

  // New Inform-based version (tjw):
  Inform fdi(NULL,0);
  setInform(fdi);

  // Scalar Field Cell-Centered------------------------------------------------
  typedef UniformCartesian<D,double> M;
  M mesh(verts);
  e_dim_tag sp[D]; for (d=0; d<D; d++) sp[d] = PARALLEL;
  CenteredFieldLayout<D,M,Cell> cl(mesh,sp,nvnodes);
  BConds<double,D,M,Cell> cbc;
  for (face=0; face < 2*D; face++) {
    if (face/2 == 0) {
      cbc[face] = new NegReflectFace<double,D,M,Cell>(face);
    } else {
      cbc[face] = new PosReflectFace<double,D,M,Cell>(face);
    }
  }
  fdi << "cbc: " << endl;
  cbc.write(fdi.getStream());
  fdi << endl;
  GuardCellSizes<D> gc(nguards);
  Field<double,D,M,Cell> A(cl,cbc,gc);
  for (d=0; d<D; d++) A[cells] += cells[d];

  fdi << endl << "--------fp3(A)-------------------------------" << endl;
  fp3(A);
  fdi << endl << "--------ggfp3(A)-----------------------------" << endl;
  ggfp3(A);
  fdi << endl << "--------agfp3(A)-----------------------------" << endl;
  agfp3(A);

  // Scalar Field Face(0)-Centered---------------------------------------------
  typedef CommonCartesianCenterings<D,1U,0U>::allFace FC;
  CenteredFieldLayout<D,M,FC> fl(mesh,sp,nvnodes);
  BConds<double,D,M,FC> fbc;
  for (face=0; face < 2*D; face++) {
    if (face/2 == 0) {
      fbc[face] = new NegReflectFace<double,D,M,FC>(face);
    } else {
      fbc[face] = new PosReflectFace<double,D,M,FC>(face);
    }
  }
  fdi << "fbc: " << endl;
  fbc.write(fdi.getStream());
  fdi << endl;
  Field<double,D,M,FC> B(fl,fbc,gc);
  NDIndex<D> faces; 
  for (d=0; d<D; d++) {
    if (d == 0) {
      faces[d] = Index(nverts[d]);
    } else {
      faces[d] = Index(nverts[d]-1);
    }
  }
  for (d=0; d<D; d++) B[faces] += faces[d];

  fdi << endl << "--------fp3(B)-------------------------------" << endl;
  fp3(B);
  fdi << endl << "--------ggfp3(B)-----------------------------" << endl;
  ggfp3(B);

  // Vektor Field vectorFace-Centered------------------------------------------
  typedef CommonCartesianCenterings<D,D>::vectorFace VFC;
  CenteredFieldLayout<D,M,VFC> vfl(mesh,sp,nvnodes);
  BConds<Vektor<double,D>,D,M,VFC> vfbc;
  for (face=0; face < 2*D; face++) {
    for (d=0; d<D; d++) {
      if (d == face/2) {
	vfbc[face*D+d] = new NegReflectFace<Vektor<double,D>,D,M,VFC>(face,d);
      } else {
	vfbc[face*D+d] = new PosReflectFace<Vektor<double,D>,D,M,VFC>(face,d);
      }
    }
  }
  fdi << "vfbc: " << endl;
  vfbc.write(fdi.getStream());
  fdi << endl;
  Field<Vektor<double,D>,D,M,VFC> C(vfl,vfbc,gc);
  // Identifiable value for unset slots (like undefined vector components in
  // vectorFace-centered Field:
  C = 9.99;

//   int component;
//   Vektor<double,D> its; for (d=0; d<D; d++) its(d) = 1.0;
//   for (component=0; component < D; component++) {
//     NDIndex<D> vfaces; 
//     for (d=0; d<D; d++) {
//       if (d == component) {
// 	vfaces[d] = Index(nverts[d]);
//       } else {
// 	vfaces[d] = Index(nverts[d]-1);
//       }
// //       assign(C[vfaces](d), vfaces[d]); // Note: op= not defined for this
//     }
//     C[vfaces] = its;
//   }
  int i,j,k;
  int counter=0;
  Vektor<double,D> vd;
  for (k=0; k < nverts[2] - 1; k++) {
    for (j=0; j < nverts[1] - 1; j++) {
      for (i=0; i < nverts[0]; i++) {
	vd = C[i][j][k].get();
	vd(0) = counter++;
	C[i][j][k] = vd;
      }
    }
  }
  counter=0;
  for (k=0; k < nverts[2] - 1; k++) {
    for (j=0; j < nverts[1]; j++) {
      for (i=0; i < nverts[0] - 1; i++) {
	vd = C[i][j][k].get();
	vd(1) = counter++;
	C[i][j][k] = vd;
      }
    }
  }
  counter=0;
  for (k=0; k < nverts[2]; k++) {
    for (j=0; j < nverts[1] - 1; j++) {
      for (i=0; i < nverts[0] - 1; i++) {
	vd = C[i][j][k].get();
	vd(2) = counter++;
	C[i][j][k] = vd;
      }
    }
  }

  setFormat(2,3);
  fdi << endl << "--------fp3(C)-------------------------------" << endl;
  fp3(C);
  fdi << endl << "--------ggfp3(C)-----------------------------" << endl;
  ggfp3(C);

  return 0;
}

/***************************************************************************
 * $RCSfile: TestFieldDebug2.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:09 $
 * IPPL_VERSION_ID: $Id: TestFieldDebug2.cpp,v 1.1.1.1 2000/11/30 21:05:09 adelmann Exp $ 
 ***************************************************************************/
