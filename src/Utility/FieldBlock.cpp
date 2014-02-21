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
 * Visit www.amas.web.psi for more details
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
#include "Utility/FieldBlock.h"
#include "Utility/IpplInfo.h"
#include "Utility/PAssert.h"
#include "Field/BrickExpression.h"
#include "Field/Field.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/LField.h"
#include "Message/Message.h"
#include "Message/Communicate.h"


#include <cstring>
#include <cstdio>

#ifdef IPPL_NETCDF
#include <netcdf.h>
#endif

//------------------------------------------------------------------
// This constructor is used to create a block of fields in a netcdf
// record. There can be an unlimited number of records. Thus, a netcdf
// file is produced which will store numFields Fields with the layout
// given into a netcdf file named fname. A new netcdf file is
// created whenever this constructor is used. 
template<class T, unsigned Dim, class Mesh, class Centering >
FieldBlock<T,Dim,Mesh,Centering>::FieldBlock(char* fname,
					     FieldLayout<Dim>& layout,
					     unsigned numFields) :
  Layout(layout), NumFields(numFields)
{
  
  strcpy(FName, fname);
  NumRecords = 0;

  // setup the file
#ifdef IPPL_NETCDF
  int Parent = 0;
  if( Ippl::Comm->myNode() == Parent ) {
    int ncid;			      // NetCDF file ID
    int rcode;

  // Create the NetCDF file, overwrite if already exists:
    ncid = nccreate(FName, NC_CLOBBER);

    // Select no-fill mode, to avoid wasting time initializing values
    // into variables with no UNLIMITED dimension upon creation:
    int fillmode = ncsetfill(ncid,NC_NOFILL); //tjwdebug
    // Initialize metadata:
    // Dimensions:
  
    const NDIndex<Dim> &domain = Layout.getDomain();
  
    // Variables:
    int dimids[4];
    if ( Dim==1 ) {
      dimids[0] = ncdimdef(ncid, "record", NC_UNLIMITED);
      dimids[1] = ncdimdef(ncid, "nx", domain[0].length());
    } 
    if ( Dim==2 ) {
      dimids[0] = ncdimdef(ncid, "record", NC_UNLIMITED);
      dimids[1] = ncdimdef(ncid, "nx", domain[0].length());
      dimids[2] = ncdimdef(ncid, "ny", domain[1].length());
    } 
    if ( Dim==3 ) {
      dimids[0] = ncdimdef(ncid, "record", NC_UNLIMITED);
      dimids[1] = ncdimdef(ncid, "nx", domain[0].length());
      dimids[2] = ncdimdef(ncid, "ny", domain[1].length());
      dimids[3] = ncdimdef(ncid, "nz", domain[2].length());
    } 
    if ( Dim>3 ) {
      ERRORMSG("FieldBlock: can't write more than 3 dimensions" << endl);
    }
    for( int i = 0 ; i < NumFields ; i++ ) {
      char varName[256];
      sprintf(varName, "var%d", i);
      ncvardef(ncid, varName, NC_DOUBLE, Dim+1, dimids);
    }
    // End (metadata) definition mode and close file (par_io requires):
    int errRet = ncendef(ncid);
    rcode = ncclose(ncid);   // Everybody/master closes file.
    if ( rcode != 0) {
      ERRORMSG("FieldBlock: ncclose() error, rcode=" << rcode << endl);
    }
  }
#endif // IPPL_NETCDF
}
//------------------------------------------------------------------
// This constructor is used to access a block of fields in a netcdf
// record. It assumes that the netcdf file fname already exists.
template<class T, unsigned Dim, class Mesh, class Centering >
FieldBlock<T,Dim,Mesh,Centering>::
FieldBlock(char* fname, FieldLayout<Dim>& layout) 
: Layout(layout)
{

    strcpy(FName, fname);
  // setup the file
#ifdef IPPL_NETCDF
  int Parent = 0;
  if( Ippl::Comm->myNode() == Parent ) {
    int i; // loop variables
    int ncid, ndims, nvars;  // NetCDF file info
    nc_type datatype;        // variable Type 
    int rcode;
    long netSize[Dim+1];
    
    // Create the NetCDF file, overwrite if already exists:
    ncid = ncopen(FName, NC_NOWRITE);
    
    // Inquire the metadata
    ncinquire(ncid, &ndims, &nvars, 0, 0);
    if( ndims != Dim+1) 
      ERRORMSG("FieldBlock: bad number of dims on read of " << FName << endl);
    NumFields = nvars;
    for( i = 0; i < Dim+1 ; i++ ) {
      ncdiminq(ncid, i, 0, &netSize[i]);
    }
    for( i = 0; i < Dim ; i++ ) {
      if( netSize[i+1] != Layout.getDomain()[i].length() ) {
	ERRORMSG("FieldBlock: encountered non-conforming size on axis ");
	ERRORMSG(i << endl);
      }
    }

    for( i = 0; i < NumFields ; i++ ) {
      ncvarinq(ncid, i, 0, &datatype, 0, 0, 0);
      if( datatype != NC_DOUBLE)
	ERRORMSG("FieldBlock: file must contain double precion data" << endl);
    }
    long numRecords;         // how many records in this file?
    ncdiminq(ncid, 0, 0, &numRecords);
    NumRecords = numRecords;
    rcode = ncclose(ncid);   // Everybody/master closes file.
    if ( rcode != 0) {
      ERRORMSG("FieldBlock: ncclose() error, rcode=" << rcode << endl);
    }
  }
#endif // IPPL_NETCDF
}
//------------------------------------------------------------------
template< class T, unsigned Dim, class Mesh, class Centering >
void FieldBlock<T,Dim,Mesh,Centering>::write(Field<T,Dim,Mesh,Centering>& f,
					     unsigned varID,
					     unsigned record)
{

#ifdef IPPL_NETCDF
  Inform msg("FieldBlock::write");
  msg.setPrintNode();

  if( varID >=NumFields ) {
    ERRORMSG(varID << " is a bad variable ID in FieldBlock::write " << endl);
    return;
  }
  int ncid;			      // NetCDF file ID
  int Parent = 0;                     // send everything to node 0

  // Open netCDF file for reading
  if( Ippl::Comm->myNode() == Parent ) {
    ncid = ncopen(FName, NC_WRITE);

    // Select no-fill mode, to avoid wasting time initializing values
    // into variables with no UNLIMITED dimension upon creation:
    int fillmode = ncsetfill(ncid,NC_NOFILL); //tjwdebug
  }
  Ippl::Comm->barrier();

  long startIndex[Dim+1];
  long countIndex[Dim+1];
  int rcode;

  int tag = Ippl::Comm->next_tag( FB_WRITE_TAG, FB_TAG_CYCLE );
  typedef LField<T,Dim>::iterator LFI;
  // ----------------------------------------
  // First loop over all the local nodes and send
  Field<T,Dim,Mesh,Centering>::iterator_if local;
  // msg << "starting to send messages: " << endl;
  for (local = f.begin_if(); local != f.end_if(); ++local) {
    // Cache some information about this local field.
    LField<T,Dim>& l = *(*local).second.get();
    NDIndex<Dim>& lo = (NDIndex<Dim>&) l.getOwned();
    NDIndex<Dim>& la = (NDIndex<Dim>&) l.getAllocated();

    // Build a message containing the owned LocalField data
    Message *mess = new Message();
    ::putMessage(*mess, lo);
    LFI msgval = l.begin();
    msgval.TryCompress();
    ::putMessage(*mess, msgval);

    // Send it.
    // msg << "sending a message from node " << Ippl::Comm->myNode();
    // msg << " to parent node " << Parent << " with tag " << tag << endl;
    Ippl::Comm->send(mess, Parent, tag);
  }
  // ----------------------------------------
  // Receive all the messages.
  // write each one to the netCDF file as it is received
  // only the Parent processor writes the data
  if( Ippl::Comm->myNode() == Parent ) {
    // we expect to receive one message from each vnode

    int numVnodes = Layout.size_iv() + Layout.size_rdv();
    // msg << " numVnodes = " << numVnodes << endl;
    
    for (int remaining = numVnodes; remaining>0; --remaining) {
      // Receive the generic message.
      int any_node = COMM_ANY_NODE;
      Message *mess = Ippl::Comm->receive_block(any_node, tag);
      PAssert(mess != 0);

      // msg << "received a message from node " << any_node;
      // msg << " on parent node " << Parent << " with tag " << tag << endl;

      // Extract the rhs BrickIterator from it.
      NDIndex<Dim> localBlock;
      LFI rhs;
      ::getMessage(*mess, localBlock);
      ::getMessage(*mess, rhs);

      // Get pointer to the data so we can free it.
      int size = 1;
      for( int i = 0 ; i < Dim ; i++ ) {
	startIndex[Dim - i] = localBlock[Dim - 1 - i].first();
	countIndex[Dim - i] = localBlock[Dim - 1 - i].length();
	size *= countIndex[Dim - i];
      }

      // unlimited dimension
      startIndex[0] = record;
      countIndex[0] = 1;
      double* buffer = new double[size]; 
      // now write the data
      int icount = 0;
      int n0,n1,n2,i0,i1,i2;
      if (rhs.IsCompressed()) {
	for (i0=0; i0<size; ++i0) buffer[icount++] = *rhs;
	// msg << " compressed value is: " << *rhs << endl;
      } else {
	switch(Dim) {
	case 1:
	  n0 = rhs.size(0);
	  for (i0=0; i0<n0; ++i0) 
	    buffer[icount++] = rhs.offset(i0);
	  break;
	case 2:
	  n0 = rhs.size(0);
	  n1 = rhs.size(1);
	  for (i0=0; i0<n0; ++i0) 
	    for (i1=0; i1<n1; ++i1) 
	      buffer[icount++] = rhs.offset(i0,i1);
	  break;
	case 3:
	  n0 = rhs.size(0);
	  n1 = rhs.size(1);
	  n2 = rhs.size(2);
	  for (i0=0; i0<n0; ++i0)
	    for (i1=0; i1<n1; ++i1)
	      for (i2=0; i2<n2; ++i2) 
		buffer[icount++] = rhs.offset(i0,i1,i2);
	  break;
	default:
	  ERRORMSG("FieldBlock: bad Dimension in Field::write()" << endl);
	  break;
	} 
      }
      // msg << " before ncvarput " << endl;
      rcode = ncvarput(ncid, varID, startIndex,countIndex, buffer);
      // msg << " after ncvarput " << endl;
      if ( rcode != 0) {
	ERRORMSG("FieldBlock: ncvarput() error, rcode=" << rcode << endl);
      }
      delete [] buffer;
      delete mess;
    }
    rcode = ncclose(ncid);   // parent closes file.
    if ( rcode != 0) {
      ERRORMSG("FieldBlock: ncclose() error, rcode=" << rcode << endl);
    }
  }
  NumRecords = NumRecords > record ? NumRecords : record;
#endif // IPPL_NETCDF
}
//--------------------------------------------------------------------
template< class T, unsigned Dim, class Mesh, class Centering >
void FieldBlock<T,Dim,Mesh,Centering>::read(Field<T,Dim,Mesh,Centering>& f,
					     unsigned varID,
					     unsigned record)
{
#ifdef IPPL_NETCDF
  Inform msg("FieldBlock::read");
  int i; // loop variables
  int ncid;                // NetCDF file info
  int Parent = 0;          // send everything to node 0

  // Create the NetCDF file, overwrite if already exists:
  if( Ippl::Comm->myNode() == Parent ) {
    ncid = ncopen(FName, NC_NOWRITE);
  }

  if( record >= NumRecords )
    ERRORMSG("invalid record on FieldBlock::read() " << endl);

  // Loop over all the Vnodes, creating an LField in each.
  long startIndex[Dim+1];
  long countIndex[Dim+1];
  int rcode;

  int tag = Ippl::Comm->next_tag( FB_READ_TAG, FB_TAG_CYCLE );
  typedef LField<T,Dim>::iterator LFI;

  // cycle through all the local vnodes and 
  // assign data to the corresponding localFields
  // only the Parent reads
  if( Ippl::Comm->myNode() == Parent ) {
    // DomainMap<NDIndex,Vnode*>::iterator for Vnodes from FieldLayout
    // this section can be reworked to straight BrickExpressions
    // the messages are being used to shake down the communications
    FieldLayout<Dim>::iterator_iv localVnode;
    for (localVnode = Layout.begin_iv() ;
	 localVnode != Layout.end_iv(); ++localVnode) {
      // Cache some information about this local vnode
      Vnode<Dim>& vn = *(*localVnode).second;
      NDIndex<Dim>& lo =  (NDIndex<Dim>&) vn.getDomain();

      int size = 1;
      for( i = 0 ; i < Dim ; i++ ) {
	startIndex[Dim - i] = lo[Dim - 1 - i].first();
	countIndex[Dim - i] = lo[Dim - 1 - i].length();
	size *= countIndex[Dim - i];
      }
      // unlimited dimension
      startIndex[0] = record;
      countIndex[0] = 1;
      double* buffer = new double[size]; 

      rcode = ncvarget(ncid, varID, startIndex, countIndex, buffer);
      if ( rcode != 0) {
	ERRORMSG("FieldBlock: ncvarget() error, rcode=" << rcode << endl);
      }
      // Build a message containing the owned LocalField data
      Message *mess = new Message();
      ::putMessage(*mess, lo);
      LFI msgval(buffer, lo, lo);
      ::putMessage(*mess, msgval);


      // Send it to the physical node
      Ippl::Comm->send(mess, vn.getNode(), tag);
      delete [] buffer;
    }
    // cycle through all the remote vnodes and 
    // send a message to the corresponding pnode
    FieldLayout<Dim>::iterator_dv remoteVnode;
    for (remoteVnode = Layout.begin_rdv() ;
	 remoteVnode != Layout.end_rdv(); ++remoteVnode) {
      // Cache some information about this remote vnode
      NDIndex<Dim>& lo =  (NDIndex<Dim>&) (*remoteVnode).first;
      Vnode<Dim>& vn        = *(*remoteVnode).second;

      int size = 1;
      for( i = 0 ; i < Dim ; i++ ) {
	startIndex[Dim - i] = lo[Dim - 1 - i].first();
	countIndex[Dim - i] = lo[Dim - 1 - i].length();
	size *= countIndex[Dim - i];
      }
      // unlimited dimension
      startIndex[0] = record;
      countIndex[0] = 1;
      double* buffer = new double[size]; 

      rcode = ncvarget(ncid, varID, startIndex,countIndex, buffer);
      if ( rcode != 0) {
	ERRORMSG("FieldBlock: ncvarget() error, rcode=" << rcode << endl);
      }

      // Build a message containing the owned LocalField data
      Message *mess = new Message();
      ::putMessage(*mess, lo);
      LFI msgval(buffer, lo, lo);
      ::putMessage(*mess, msgval);

      // Send it to the physica node
      Ippl::Comm->send(mess, vn.getNode(), tag);
      delete [] buffer;
    }
    rcode = ncclose(ncid);   // Parent closes file.
    if ( rcode != 0) {
      ERRORMSG("FieldBlock: ncclose() error, rcode=" << rcode << endl);
    }
  }
  
  // receive the messages
  int numLocalVnodes = Layout.size_iv();
  for (int remaining = numLocalVnodes; remaining>0; --remaining) {
    // Receive the generic message.
    int any_node = COMM_ANY_NODE;
    Message *mess = Ippl::Comm->receive_block(any_node, tag);
    PAssert(mess != 0);

    // Extract the rhs BrickIterator from it.
    NDIndex<Dim> localBlock;
    LFI rhs;
    ::getMessage(*mess, localBlock);
    ::getMessage(*mess, rhs);

    // to which local vnode does this correspond?
    LField<T,Dim>* myLField;
    bool flag = 1;
    // DomainMap<NDIndex,LField*>::iterator for LocalFields
    Field<T,Dim,Mesh,Centering>::iterator_if local;
    for (local = f.begin_if(); local != f.end_if(); ++local) {
      myLField = (*local).second.get();
      const NDIndex<Dim>& lo = myLField->getOwned();
      if( lo == localBlock ) {
	flag = 0;
	break; 
      }
    }
    if(flag) 
      ERRORMSG("FieldBlock::read: did not match the local NDIndex" << endl);
    // now read the data into the LocalFields
    // Build the lhs brick iterator.
    LFI lhs = myLField->begin();
    int n0,n1,n2,i0,i1,i2;
    if( lhs.IsCompressed() && rhs.IsCompressed() ) {
      *lhs = *rhs;
    }
    if( !lhs.IsCompressed() && rhs.IsCompressed() ) {
      for (i0=0; i0 < localBlock.size(); ++i0) {
	*lhs = *rhs;
	++lhs;
      }
    }
    if( !rhs.IsCompressed() ) {
      if( lhs.IsCompressed() ) {
	myLField->Uncompress();
	lhs = myLField->begin();
      }
      switch(Dim) {
      case 1:
	n0 = rhs.size(0);
	for (i0=0; i0<n0; ++i0) {
	  lhs.offset(i0) = *rhs;
	  ++rhs;
	}
	break;
      case 2:
	n0 = rhs.size(0);
	n1 = rhs.size(1);
	for (i0=0; i0<n0; ++i0) 
	  for (i1=0; i1<n1; ++i1) {
	    lhs.offset(i0,i1) = *rhs;
	    ++rhs;
	  }
	break;
      case 3:
	n0 = rhs.size(0);
	n1 = rhs.size(1);
	n2 = rhs.size(2);
	for (i0=0; i0<n0; ++i0) 
	  for (i1=0; i1<n1; ++i1)
	    for (i2=0; i2<n2; ++i2) {
	      lhs.offset(i0,i1,i2) = *rhs;
	      ++rhs;
	    }
	break;
      default:
	ERRORMSG("FieldBlock: bad Dimension in write()" << endl);
	break;
      }
    }
    // Do the assignment.
    //    BrickExpression<Dim,BI,BI,AppAssign<T,T> > (lhs,rhs).apply();
    // Free the memory.
    delete mess;
  } // loop over the vnodes

#endif // IPPL_NETCDF
}
//--------------------------------------------------------------------
/***************************************************************************
 * $RCSfile: FieldBlock.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: FieldBlock.cpp,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
