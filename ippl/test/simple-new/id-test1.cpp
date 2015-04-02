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
 ***************************************************************************/

/***************************************************************************
 * id-test1.cpp - a simple test program that
 */

#include "Ippl.h"

template<class T, unsigned Dim>
class PDBAtoms : public ParticleBase< ParticleInteractLayout<T, Dim> > {

typedef  ParticleInteractLayout<T, Dim> Layout_t;

public:
  // the constructor: in this case, this class will determine the field
  // layout from the particle data, and so the only argument is the name
  // of the PDB file.

  PDBAtoms() {
      // initialize the base class, by creating a new layout
      initialize(new Layout_t());
      this->create(10);
      this->update();
  }
    
};



// the main routine
int main(int argc, char *argv[]) {

  // initialize Ippl
  Ippl ippl(argc, argv);
  Inform PDBmsg(argv[0], INFORM_ALL_NODES);

  PDBAtoms<float,3> *pdb = new PDBAtoms<float,3>();

  PDBmsg << "Reading complete.  Total atoms = " << pdb->getTotalNum();
  PDBmsg << ", local = " << pdb->getLocalNum() << endl;

  for (int i=0; i<pdb->getLocalNum(); i++)
        PDBmsg << "ID = " << pdb->ID[i] << endl;

  PDBmsg << endl;

  pdb->destroy(5,2);
  pdb->update();

  for (int i=0; i<pdb->getLocalNum(); i++)
        PDBmsg << "ID = " << pdb->ID[i] << endl;

  pdb->createWithID(5);
  PDBmsg << "create doene ... = " << endl;
  pdb->update();


  PDBmsg << endl;
  for (int i=0; i<pdb->getLocalNum(); i++)
        PDBmsg << "ID = " << pdb->ID[i] << endl;


  Ippl::Comm->barrier();	// wait for all nodes to finish
  return 0;
}

/***************************************************************************
 * $RCSfile: pdbtest.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: pdbtest.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
