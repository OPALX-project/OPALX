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

#include "Ippl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string>
#include <fstream>

const unsigned Dim = 3;

typedef Cartesian<Dim,double> MT; //typedef Vektor<double,Dim> T;
typedef double T;


 

// dimension of our positions

void dumpVTK(Field<T,Dim,MT,Vert> &EFD, NDIndex<Dim> lDom, int nx, int ny, int nz, int iteration, double dx, double dy, double dz) {

    ofstream vtkout;
    vtkout.precision(10);
    vtkout.setf(ios::scientific, ios::floatfield);

    std::stringstream fname;
    fname << "data/A_";
    fname << setw(4) << setfill('0') << iteration;
    fname << ".vtk";

    //SERIAL at the moment
    //if (Ippl::myNode() == 0) {

    // open a new data file for this iteration
    // and start with header
    vtkout.open(fname.str().c_str(), ios::out);
    vtkout << "# vtk DataFile Version 2.0" << endl;
    vtkout << "cref" << endl;
    vtkout << "ASCII" << endl;
    vtkout << "DATASET STRUCTURED_POINTS" << endl;
    vtkout << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    vtkout << "ORIGIN 0.0 0.0 0.0" << endl;
    vtkout << "SPACING " << dx << " " << dy << " " << dz << endl;
    vtkout << "POINT_DATA " << nx*ny*nz << endl;

    vtkout << "VECTORS A-Field float" << endl;
    for(int z=lDom[2].first(); z<lDom[2].last(); z++) {
      for(int y=lDom[1].first(); y<lDom[1].last(); y++) {
	for(int x=lDom[0].first(); x<lDom[0].last(); x++) {
	  T tmp = EFD[x][y][z].get();
	  vtkout << tmp << "\t" << tmp << "\t" << tmp << endl;
	}
      }
    }
    // close the output file for this iteration:
    vtkout.close();
}


int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  Inform msg2all(argv[0],INFORM_ALL_NODES);

  msg << "Core Refinement Test" << endl;

  int i, j;
    
  /*
    Construct the Domain (0,0) - (1,1)
    with a fine domain from 0.2 - 0.3

  */

  double meshspacingC = 0.1;
  double meshspacingF = meshspacingC / 10.;

  int nx = (int) (1.0 / meshspacingC) + (0.1 /  meshspacingF);
  int ny = (int) (1.0 / meshspacingC) + (0.1 /  meshspacingF);
  int nz = (int) (1.0 / meshspacingC) + (0.1 /  meshspacingF);

  msg << "nx=ny= " << nx << endl;


  // Sizes:
  unsigned nverts[Dim], ncells[Dim];
  unsigned totverts=1, totcells=1;
  int d;
  for (d=0; d<Dim; d++) {
      ncells[d] = nx - 1;
      nverts[d] = nx;
      totcells *= ncells[d];
      totverts *= nverts[d];
  }
  
  NDIndex<Dim> verts, cells;
  for (d=0; d<Dim; d++) {
      verts[d] = Index(nverts[d]);
      cells[d] = Index(ncells[d]);
  }
  
  double* h[Dim];

  for (d=0; d<Dim; d++) 
      h[d] = new double[nverts[d]];
  

  Vektor<double,Dim> origin;
  for (int d=0; d<Dim; d++) 
      origin(d) = 0.0; //d + 1.0;
  
  double hx = 1.0/nx;
  for (int d=0; d<Dim; d++) {
      for (int vert=0; vert < nverts[d]; vert++) {
	  (h[d])[vert] = hx;
      }
  }
  // Mesh boundary conditions:
  MeshBC_E mbc[2*Dim];
  for (unsigned b=0; b < (2*Dim); b++) 
      mbc[b] = Reflective;

  // Test constructing mesh, and then setting spacing, origin, BC's


  MT mesh(verts);
  mesh.set_origin(origin);
  mesh.set_meshSpacing(h);
  mesh.set_MeshBC(mbc);
  mesh.storeSpacingFields();


  // Construct CenteredFieldLayout's using this for Vert and Cell centering:
  e_dim_tag edt[Dim];
  for (int d=0; d<Dim; d++) 
      edt[d] = PARALLEL;
  CenteredFieldLayout<Dim,MT,Cell> cFL(mesh, edt);
  CenteredFieldLayout<Dim,MT,Vert> vFL(mesh, edt);

  // Use 1 guard layer in all Field's:
  GuardCellSizes<Dim> gc(1);

  // Use linear negative reflecting conditions:

  // Vectors:
  BConds<T,Dim,MT,Vert> vvbc;
  BConds<T,Dim,MT,Cell> vcbc;
  // Scalars:
  BConds<double,Dim,MT,Cell> scbc; //  BConds<double,D,M,Vert> svbc;

  for (int face=0; face<2*Dim; face++) {
      vvbc[face]  = new NegReflectFace<T,Dim,MT,Vert>(face);
      vcbc[face]  = new NegReflectFace<T,Dim,MT,Cell>(face);
      //    svbc[d] = new NegReflectFace<double,D,M,Vert>(face);
      scbc[face]  = new NegReflectFace<double,Dim,MT,Cell>(face);
  }

  Field<T,Dim,MT,Vert> A(mesh, vFL, gc, vvbc); // aVectorFieldVert(mesh, vFL, gc, vvbc);
  Field<T,Dim,MT,Vert> B(mesh, vFL, gc, vvbc); // aVectorFieldVert(mesh, vFL, gc, vvbc);
  NDIndex<Dim> lDom = vFL.getLocalNDIndex();

  int iterations = atoi(argv[1]);
  
  A = 0.0;
  B = 0.0;
  
  int centerX = nx/2;
  int centerY = ny/2;
  int centerZ = nz/2;
  
  A[centerX][centerY][centerZ] = T(1.0*iterations);

  Index I(cells[0]);
  Index J(cells[1]);
  Index K(cells[2]);

  double fact = 1.0/9.0;

  for(int iter = 0 ; iter < iterations ; iter++ ) {
      
    B[I][J][K]  = fact * (
			  A[I  ][J  ][K+1] + A[I  ][J  ][K-1] +
			  A[I  ][J+1][K  ] + A[I  ][J-1][K  ] +
			  A[I+1][J  ][K  ] + A[I-1][J  ][K  ]
			  );
    B[I][J][K] += fact * (
			  A[I+1][J+1][K  ] + A[I+1][J-1][K  ] +
			  A[I  ][J  ][K  ] + A[I-1][J+1][K  ] +
			  A[I-1][J-1][K  ]
			  );
    B[I][J][K] += fact * (
			  A[I+1][J+1][K+1] + A[I+1][J  ][K+1] +
			  A[I+1][J-1][K+1] + A[I  ][J+1][K+1] +
			  A[I  ][J-1][K+1] + A[I-1][J+1][K+1] +
			  A[I-1][J  ][K+1] + A[I-1][J-1][K+1]
                  );
    B[I][J][K] += fact * (
			  A[I+1][J+1][K-1] + A[I+1][J  ][K-1] +
			  A[I+1][J-1][K-1] + A[I  ][J+1][K-1] +
			  A[I  ][J-1][K-1] + A[I-1][J+1][K-1] +
			  A[I-1][J  ][K-1] + A[I-1][J-1][K-1]
			  );
    A = B;
    
    dumpVTK(A,lDom,nx,ny,nz,iter,hx,hx,hx);
    
    cout << "iter = " << iter << " sum = " << sum(A) << endl;
      
  }
  msg << "Core Refinement Test End." << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: spatial.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: spatial.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/


/*
 NDIndex<Dim> lDom = vFL.getLocalNDIndex();

  Index ii = lDom[0];
  Index jj = lDom[1];

  for (int i = ii.min(); i < ii.max(); i++) {
      for (int j = jj.min(); j < jj.max(); j++) {
	  NDIndex<Dim> el(Index(i,i),Index(j,j));
	 msg << i << "  " << j << "   " << aField.localElement(el) << " V= " << mesh.getCellVolume(el) << " dh= " << mesh.getDeltaVertex(el) << endl;
      }
  }

  msg << aField << endl;
  msg << mesh << endl;

  msg << "Core Refinement Test End." << endl;
  return 0;

*/
