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

/***************************************************************************
 * pdbtest.cpp - a simple test program that reads in a PDB file, initializes
 * a user particle object, and performs test operations with this data.
 *
 * Note: the PDB file format is a standard method for storing biopolymer atomic
 * coordinate data (such as crystal structures of proteins and nucleic acids).
 * See http://www.pdb.bnl.gov/ for more info on PDB files.
 ***************************************************************************/

#include "Ippl.h"
#include <stdio.h>

//
// user-defined class which will store a set of atoms read from a PDB file
//
template<class T, unsigned Dim>
class PDBAtoms : public ParticleBase< ParticleInteractLayout<T, Dim> > {

public:
  // simple class to store an atom name
  struct AtomName {
    char name[8];
    AtomName() { name[0] = '\0'; }
    void setName(char *newname) { strcpy(name, newname); }
    Message& putMessage(Message& m) {
      ::putMessage(m, name, name + 8);
      return m;
    }
    Message& getMessage(Message& m) {
      ::getMessage(m, name);
      return m;
    }
  };

public:
  // the attributes for this set of particles (atoms).
  ParticleInteractAttrib<T> Beta;
  ParticleInteractAttrib<T> Occup;
  ParticleAttrib<AtomName> Name;
  ParticleAttrib<AtomName> ResName;
  ParticleAttrib<int> ResID;
  ParticleInteractAttrib<unsigned> Num;

  // the constructor: in this case, this class will determine the field
  // layout from the particle data, and so the only argument is the name
  // of the PDB file.
  PDBAtoms(char *fname, double interrad) {
    // initialize the base class, by creating a new layout
    initialize(new Layout_t());

    // register our attributes with the base class
    addAttribute(Beta);
    addAttribute(Occup);
    addAttribute(Name);
    addAttribute(ResName);
    addAttribute(ResID);
    addAttribute(Num);

    // set the interaction radius for the atoms.
    getLayout().setInteractionRadius(interrad);

    // read in the file; this will also do an update to properly disburse
    // the atoms among the processors
    read_pdb_file(fname);
  }

private:
  // enumeration with types of PDB records
  enum { PDB_UNKNOWN, PDB_ATOM, PDB_REMARK, PDB_END };

  // The routine which reads the PDB file and adds the data to our storage.
  // This will only store every Nth atom on the local node, where N is the
  // number of processors.  After all atoms are read in, 'update'
  // is invoked to make sure all the particles are moved to their proper
  // node (based on their position and the FieldLayout).  update will also
  // calculate the spatial layout for the particles, based on the maximum
  // extent of all the particles' locations.
  void read_pdb_file(char *fname) {

    // find out our local processor ID number
    int myN  = Ippl::myNode();	
    int totN = Ippl::getNodes();

    // open the file; if we can't, just quit now
    FILE *f = fopen(fname, "r");
    if (f == 0) {
      ERRORMSG("Cannot open PDB file " << fname << endl);
      return;
    }

    // scan through the file, adding every myNth atom
    unsigned totalatoms = 0;
    unsigned localatoms = 0;
    char record[128];		// PDB records are never longer than 80 chars
    char name[8], resname[8], resid[8];
    float pos[3], occup, beta;
    int rectype, anum;
    while ((rectype = read_pdb_record(f, record)) != PDB_END) {
      if (rectype == PDB_ATOM) {
	// add it to our local atoms if we need to
	if ((totalatoms++ % totN) == myN) {
	  // break an atom record into separate fields
	  anum = get_pdb_fields(record, name, resname, resid,
				pos[0], pos[1], pos[2], occup, beta);
	  create(1);
	  R[localatoms] = Vektor<T,Dim>(pos[0], pos[1], pos[2]);
	  Beta[localatoms] = beta;
	  Occup[localatoms] = occup;
	  Num[localatoms] = anum;
	  ResID[localatoms] = atoi(resid);
	  Name[localatoms].setName(name);
	  ResName[localatoms].setName(resname);
	  localatoms++;
	}

	// we've finished processing this atom
	totalatoms++;
      }
    }

    // close the PDB file now
    fclose(f);

    // finally, do an update to move particles around
    update();
  }

  // read the next record from the PDB file, and return the type of record
  // found.
  int read_pdb_record(FILE *f, char *retStr) {
    char inbuf[128], inbuf2[128];
    char *tok;
  
    // read the next line
    if(inbuf != fgets(inbuf, 127, f)) {
      strcpy(retStr,"");
      return PDB_END;
    } else {
      // initially, we return an empty string
      strcpy(retStr,"");

      // remove the newline character, if there is one
      if(inbuf[strlen(inbuf)-1] == '\n')
	inbuf[strlen(inbuf)-1] = '\0';
      strcpy(inbuf2,inbuf);
      tok = strtok(inbuf2, " \t");
      if(tok == NULL) {
	return PDB_UNKNOWN;
      } else if(strcmp(tok,"REMARK") == 0) {
	return PDB_REMARK;
      } else if(strcmp(tok,"ATOM") == 0 || strcmp(tok,"HETATM") == 0) {
	strcpy(retStr, inbuf);
	return PDB_ATOM;
      } else if(strcmp(tok,"END") == 0) {
	return PDB_END;
      } else {
	return PDB_UNKNOWN;
      }
    }
  }

  // Break a pdb ATOM record into it's fields.  The user must provide the
  // necessary space to store the atom name, residue name, and segment name.
  // Character strings will be null-terminated. Returns the atom serial number.
  int get_pdb_fields(char *record, char *name, char *resname, char *resid,
		     float& x, float& y, float& z, float& occup, float& beta) {
    int i, len, num;
    char numstr[9];
    numstr[8] = '\0';

    // get serial number
    sscanf(record + 6, "%d", &num);

    // get atom name
    strncpy(name, record + 12, 4);
    name[4] = '\0';
    while((len = strlen(name)) > 0 && name[len-1] == ' ')
      name[len-1] = '\0';
    while(len > 0 && name[0] == ' ') {
      for(i=0; i < len; i++)  name[i] = name[i+1];
      len--;
    }

    // get residue name
    strncpy(resname, record + 17, 3);
    resname[3] = '\0';
    while((len = strlen(resname)) > 0 && resname[len-1] == ' ')
      resname[len-1] = '\0';
    while(len > 0 && resname[0] == ' ') {
      for(i=0; i < len; i++)  resname[i] = resname[i+1];
      len--;
    }

    // get residue id number
    strncpy(resid, record + 22, 4);
    resid[4] = '\0';
    while((len = strlen(resid)) > 0 && resid[len-1] == ' ')
      resid[len-1] = '\0';
    while(len > 0 && resid[0] == ' ') {
      for(i=0; i < len; i++)  resid[i] = resid[i+1];
      len--;
    }

    // get coordintes, occupancy, and beta-factor
    strncpy(numstr, record + 30, 8);
    x = atof(numstr);
    strncpy(numstr, record + 38, 8);
    y = atof(numstr);
    strncpy(numstr, record + 46, 8);
    z = atof(numstr);
    strncpy(numstr, record + 54, 6);
    occup = atof(numstr);
    strncpy(numstr, record + 60, 6);
    beta = atof(numstr);

    // return atom serial number
    return num;
  }

};



// the main routine
int main(int argc, char *argv[]) {

  // initialize Ippl
  Ippl ippl(argc, argv);
  Inform PDBmsg(argv[0], INFORM_ALL_NODES);

  // get PDB filename from the command line
  if (argc < 2 || argc > 3) {
    ERRORMSG("Usage: " << argv[0] << " <pdb filename> [radius]" << endl);
    exit(1);
  }

  char *fname = argv[1];
  float radius = (argc > 2 ? atof(argv[2]) : 3.5);
  if (radius <= 0.0)
    radius = 3.5;

  // Create storage for PDB atoms, by creating a new PDBAtom instance.
  // The PDBAtoms constructor will read the data file and initialize the
  // particle data on all the nodes (including a call to the update()
  // routine).  The second argument to the PDBAtoms constructor is the
  // interaction radius for the particles (here, the same value is used
  // for all particles).
  PDBmsg << "Reading coordinates from PDB file " << argv[1] << " ..." << endl;
  PDBAtoms<float,3> pdb(fname, radius);
  PDBmsg << "Reading complete.  Total atoms = " << pdb.getTotalNum();
  PDBmsg << ", local = " << pdb.getLocalNum() << endl;

  // add periodic BC's to the particles
  pdb.getBConds()[0] = ParticlePeriodicBCond;
  pdb.getBConds()[1] = ParticlePeriodicBCond;
  pdb.getBConds()[2] = ParticlePeriodicBCond;
  pdb.getBConds()[3] = ParticlePeriodicBCond;
  pdb.getBConds()[4] = ParticlePeriodicBCond;
  pdb.getBConds()[5] = ParticlePeriodicBCond;

  // loop through the particles, and find the nearest neighbors of each one
  for (unsigned i=0; i < pdb.getLocalNum(); ++i) {
    // create empty iterators used to access the nearest neighbors
    PDBAtoms<float,3>::pair_iterator plisti;
    PDBAtoms<float,3>::pair_iterator plistf;

    // fill these iterators with data for the ith particle's neighbors.
    // after this call, plisti points to the start of an array of the local
    // indices of the neighbors of particle i, while plistf points just
    // past the end of this array.  If plisti == plistf, this particle
    // has no neighbors.
    pdb.getLayout().getPairlist(i, plisti, plistf, pdb);

    PDBmsg << "For ID = " << pdb.ID[i] << ", num neighbors = ";
    PDBmsg << (plistf - plisti) << ", R = " << pdb.R[i] << endl;

    // perform calculation for each neighbor of this particle.
    for ( ; plisti != plistf ; ++plisti ) {

      // get local index and sep^2 of neighbor
      PDBAtoms<float,3>::pair_t nbrData = *plisti;

      // print out their values
      PDBmsg << "  Neighbor:";
      PDBmsg <<  " ID="             << pdb.ID[nbrData.first];
      PDBmsg << ", sep^2="          << nbrData.second;
      PDBmsg << ", pos="            << pdb.R[nbrData.first] << endl;
    }
  }

  Ippl::Comm->barrier();	// wait for all nodes to finish
  return 0;
}

/***************************************************************************
 * $RCSfile: pdbtest.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:03:43 $
 * IPPL_VERSION_ID: $Id: pdbtest.cpp,v 1.1.1.1 2000/11/30 21:03:43 adelmann Exp $ 
 ***************************************************************************/
