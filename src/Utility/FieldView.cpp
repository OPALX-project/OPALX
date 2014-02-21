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
#include "Utility/FieldView.h"
#include "Utility/IpplInfo.h"
#include "Utility/PAssert.h"
#include "Field/BrickExpression.h"
#include "Field/Field.h"
#include "Field/LField.h"
#include "Message/Message.h"


#ifdef IPPL_GDL
#include <gdl.h>
#endif

//------------------------------------------------------------------
// attach a 2D Field to a FieldView
template<class T, unsigned Dim, class Mesh, class Centering >
FieldView<T,Dim,Mesh,Centering>::FieldView(Field<T,Dim,Mesh,Centering>& f, 
					   unsigned scaleX, 
					   unsigned scaleY, 
					   unsigned minSizeX,
					   unsigned minSizeY,
					   unsigned parent) : 
  MyField(f), SliceDim(0),ScaleX(scaleX), ScaleY(scaleY), 
  MinSizeX(minSizeX), MinSizeY(minSizeY), Parent(parent) 
{ 

#ifdef IPPL_GDL
  if (Ippl::Comm->myNode() == Parent) {
    dummy = GDL_OpenDisplay(X11FB);
    GDL_SetColormap(RAINBOW_BLUE,NULL);
    PInsist(Dim == 2,
            "This FieldView constructor only works for a 2D Field!!");
    // construct the LField which will be used to coalesce the data
    NDIndex<2U> sliceDomain;
    sliceDomain[0] = MyField.getDomain()[0];
    sliceDomain[1] = MyField.getDomain()[1];
    MyLField = new LField<T,2U>(sliceDomain, sliceDomain);
    // since all LFields are born compressed, we need
    // to explicitly uncompress the field
    MyLField->Uncompress();
    init_map();
  }
#endif // IPPL_GDL
}
//------------------------------------------------------------------
// attach a 3D Field to a FieldView 
template<class T, unsigned Dim, class Mesh, class Centering >
FieldView<T,Dim,Mesh,Centering>::FieldView(unsigned sliceDim,
					   Field<T,Dim,Mesh,Centering>& f, 
					   unsigned scaleX, 
					   unsigned scaleY, 
					   unsigned minSizeX,
					   unsigned minSizeY,
					   unsigned parent) : 
  MyField(f), SliceDim(sliceDim),ScaleX(scaleX), ScaleY(scaleY), 
  MinSizeX(minSizeX), MinSizeY(minSizeY), Parent(parent) 
{

#ifdef IPPL_GDL
  if (Ippl::Comm->myNode() == Parent) {
    dummy = GDL_OpenDisplay(X11FB);
    GDL_SetColormap(RAINBOW_BLUE,NULL);
    PInsist(SliceDim < Dim,
            "Invalid slice dimension specified in FieldView constructor!!"); 
    PInsist(Dim == 3,
            "This FieldView constructor only works for a 3D Field!!");
    // construct the LField which will be used to coalesce the data
    NDIndex<2U> sliceDomain;
    unsigned ix = SliceDim < 1 ? 1 : 0;
    unsigned iy = SliceDim < 2 ? 2 : 1;
    sliceDomain[0] = MyField.getDomain()[ix];
    sliceDomain[1] = MyField.getDomain()[iy];
    MyLField = new LField<T,2U>(sliceDomain, sliceDomain);
    // since all LFields are born compressed, we need
    // to explicitly uncompress the field
    MyLField->Uncompress();
    init_map();
  }
#endif // IPPL_GDL
}
//------------------------------------------------------------------
template<class T, unsigned Dim, class Mesh, class Centering >
FieldView<T,Dim,Mesh,Centering>::~FieldView() 
{
  
  
#ifdef IPPL_GDL
  delete MyLField; 
  delete [] MapX;
  delete [] MapY;
  delete [] Data;
#endif // IPPL_GDL
}
//------------------------------------------------------------------
// view the 2D data in a GDL window
template<class T, unsigned Dim, class Mesh, class Centering >
void FieldView<T,Dim,Mesh,Centering>::void_view(int &r) 
{
  
  
#ifdef IPPL_GDL
  update_2D_data();
  /* only 0 should do apply map */
  if (Ippl::Comm->myNode() == Parent) 
    r = apply_map();
  else
    r = 1;
#endif // IPPL_GDL
}
//------------------------------------------------------------------
// view a slice of a 3D field in a GDL window
template<class T, unsigned Dim, class Mesh, class Centering >
void FieldView<T,Dim,Mesh,Centering>::void_view(unsigned slice, int &r) 
{
  
  

#ifdef IPPL_GDL
  update_3D_data(slice);
  /* only 0 should do apply map */
  if (Ippl::Comm->myNode() == Parent)
    r = apply_map();
  else
    r = 1;
#endif // IPPL_GDL
}
//------------------------------------------------------------------
// view a slice of a 3D field in a GDL window
template<class T, unsigned Dim, class Mesh, class Centering >
void FieldView<T,Dim,Mesh,Centering>::void_apply_map(int &r)
{

#ifdef IPPL_GDL
  int icount = 0;
  LField<T,2U>::iterator liter = MyLField->begin();
  // make the origin at the lower left
  int domainY = MyField.getLayout().getDomain()[1].length();
  for( int j = 0 ; j < SizeY ; j++) {
    for( int i = 0 ; i < SizeX ; i++) {
      Data[icount++] = liter.offset(MapX[i],domainY-MapY[j]-1);
    }
  }
#endif // IPPL_GDL
#ifdef IPPL_GDL
  r = GDL_Display_double(Data, SizeX, SizeY);
#else
  r = dummy;
#endif // IPPL_GDL
}
//------------------------------------------------------------------
// draw all the data together onto the Parent process for viewing
template<class T, unsigned Dim, class Mesh, class Centering >
void 
FieldView<T,Dim,Mesh,Centering>::update_2D_data(void) 
{

#ifdef IPPL_GDL
  int tag = Ippl::Comm->next_tag( FV_2D_TAG, FV_TAG_CYCLE );
  typedef LField<T,Dim>::iterator LFI;
  
  // ----------------------------------------
  // First loop over all the local nodes and send
  Field<T,Dim,Mesh,Centering>::iterator_if local;
  for (local = MyField.begin_if(); local != MyField.end_if(); ++local) {
    // Cache some information about this local field.
    LField<T,Dim> &l = *(*local).second;
    NDIndex<Dim>& lo = (NDIndex<Dim>&) l.getOwned();
    NDIndex<Dim>& la = (NDIndex<Dim>&) l.getAllocated();
    l.Uncompress();
    T* lp = l.getP();

    // Build a message containing the owned LocalField data
    Message *mess = new Message();
    ::putMessage(*mess, lo);
    LFI msgval(lp,lo,la);
    ::putMessage(*mess, msgval);

    // Send it.
    Ippl::Comm->send(mess, Parent, tag);
  }
  
  // ----------------------------------------
  // Receive all the messages.
  if( Ippl::Comm->myNode() == Parent ) {
    // we expect to receive one message from each vnode
    int numVnodes = MyField.getLayout().size_iv() + 
      MyField.getLayout().size_rdv();
    for (int remaining = numVnodes; remaining>0; --remaining) {
      // Receive the generic message.
      int any_node = COMM_ANY_NODE;
      Message *mess = Ippl::Comm->receive_block(any_node, tag);
      PAssert(mess != 0);

      // Extract the rhs BrickIterator from it.
      NDIndex<Dim> localBlock;
      T rhs_compressed_data;
      LFI rhs(rhs_compressed_data);
      localBlock.getMessage(*mess);
      rhs.getMessage(*mess);

      // Build the lhs brick iterator.
      LFI lhs(MyLField->getP(), localBlock, MyLField->getAllocated());

      // Do the assignment.
      BrickExpression<Dim,LFI,LFI,OpAssign > (lhs,rhs).apply();

      // Free the memory.
      delete mess;
    }
  }
#endif // IPPL_GDL
}
//------------------------------------------------------------------
template<class T, unsigned Dim, class Mesh, class Centering >
void 
FieldView<T,Dim,Mesh,Centering>::update_3D_data(unsigned slice) 
{

#ifdef IPPL_GDL
  Inform testmsg("FieldView");
  PInsist(Dim == 3,
          "FieldView::update_3D_data is only valid for a 3D Field!!");
  int tag = Ippl::Comm->next_tag( FV_3D_TAG, FV_TAG_CYCLE );
  typedef LField<T,3U>::iterator LFI;
  
  const Index& sl = MyField.getDomain()[SliceDim];
  if ((slice < sl.first()) || (slice >= sl.first() + sl.length())) {
    ERRORMSG("FieldView: bad slice choice of "<< slice << " in dim ");
    ERRORMSG(SliceDim << endl);
    return;
  }
  // find the domain which represents the 2D slice plane
  NDIndex<3U> sliceDomain(MyField.getDomain());
  sliceDomain[SliceDim] = Index(slice,slice);
  //  testmsg << " Slice is " << sliceDomain << endl;
  // ----------------------------------------
  // First loop over all the local nodes and send
  Field<T,3U,Mesh,Centering>::iterator_if local;
  int vnode = 0;
  for (local = MyField.begin_if(); local != MyField.end_if(); ++local) {
    //    testmsg << " vnode = " << vnode << endl;
    // Cache some information about this local field.
    // testmsg << " vnode is " << vnode << endl;
    LField<T,3U> &l = *(*local).second;
    NDIndex<3U>& lo = (NDIndex<3U>&) l.getOwned();
    NDIndex<3U>& la = (NDIndex<3U>&) l.getAllocated();
    l.Uncompress();
    T* lp = l.getP();

    // find the intersection with the slice
    NDIndex<3U> intersection = lo.intersect( sliceDomain );
    // testmsg << " intersection is " << intersection << endl;

    // Build a message containing a slice of the owned LocalField data
    Message *mess = new Message();
    //    mess->put(intersection);
    ::putMessage(*mess, intersection);
    // testmsg << " sending: " << endl;
    // testmsg << " lo = " << intersection << endl;
    // testmsg << " la = " << la << endl;
    LFI msgval(lp,intersection,la);

    //    LFI data = l.begin();
    // print out the whole field:
    //    int i,j,k;
    //    for( i = 0 ; i < data.size(0) ; i++) {
    //      for( j = 0 ; j < data.size(1) ; j++) {
    //	for( k = 0 ; k < data.size(2) ; k++) {
    //	  testmsg << "data["<<i<<"]["<<j<<"]["<<k<<"]= "<<data.offset(i,j,k)<<endl;
    //	}
    //      }
    //    }

    //    mess->put(msgval);
    ::putMessage(*mess, msgval);

    //    vnode++;
    // Send it.
    Ippl::Comm->send(mess, Parent, tag);
  }
  
  // ----------------------------------------
  // Receive all the messages.
  if( Ippl::Comm->myNode() == Parent ) {
    // we expect to receive one message from each vnode
    int numVnodes = MyField.getLayout().size_iv() + 
      MyField.getLayout().size_rdv();
    for (int remaining = numVnodes; remaining>0; --remaining) {
      // Receive the generic message.
      int any_node = COMM_ANY_NODE;
      Message *mess = Ippl::Comm->receive_block(any_node, tag);
      PAssert(mess != 0);

      // Extract the rhs BrickIterator from it.
      NDIndex<3U> localBlock;
      T rhs_compressed_data;
      LFI rhs(rhs_compressed_data);
      localBlock.getMessage(*mess);
      rhs.getMessage(*mess);

      unsigned ix = SliceDim < 1 ? 1 : 0;
      unsigned iy = SliceDim < 2 ? 2 : 1;

      LField<double,2U>::iterator liter = MyLField->begin();

      if(rhs.size(SliceDim) > 0) {
	int f0 = localBlock[ix].first();
	int f1 = localBlock[iy].first();
	int n0 = localBlock[ix].length();
	int n1 = localBlock[iy].length();
	//	cout << " f0 = " << f0;
	//	cout << " f1 = " << f1;
	//	cout << " n0 = " << n0;
	//	cout << " n1 = " << n1;

	for (int i1=0; i1<n1; ++i1)
	  for (int i0=0; i0<n0; ++i0) {
	    //	    cout << "rhs["<<i0<<"]["<<i1<<"]= "<< *rhs << endl;
	      liter.offset(f0+i0,f1+i1) = *rhs;
	      ++rhs;
	    }
      }

      // Free the memory.
      delete mess;
    }
  }
#endif // IPPL_GDL
}
//------------------------------------------------------------------
// helper functions to fit the data into the viewing port
template<class T, unsigned Dim, class Mesh, class Centering >
void 
FieldView<T,Dim,Mesh,Centering>::init_map(void) 
{
  
  
#ifdef IPPL_GDL
  int i,j;

  int domainX, domainY;
  unsigned ix, iy;
  switch(Dim) {
  case 2:
    domainX = MyField.getLayout().getDomain()[0].length();
    domainY = MyField.getLayout().getDomain()[1].length();
    break;
  case 3:
    ix = SliceDim < 1 ? 1 : 0;
    iy = SliceDim < 2 ? 2 : 1;
    domainX = MyField.getLayout().getDomain()[ix].length();
    domainY = MyField.getLayout().getDomain()[iy].length();
    break;
  default:
    ERRORMSG("FieldView: bad dimension " << Dim << " for FieldView mapping");
    ERRORMSG(endl);
    break;
  }
  // set scale factors in x and y directions
  SizeX = domainX;
  while(SizeX < MinSizeX) {
    ScaleX++;
    SizeX = ScaleX * domainX;
  }
  SizeY = domainY;
  while(SizeY < MinSizeY) {
    ScaleY++;
    SizeY = ScaleY * domainY;
  }
  // allocate memory
  MapX = new int[SizeX];
  MapY = new int[SizeY];
  Data = new T[SizeX*SizeY];
  
  int icount = 0;
  for( i = 0 ; i < domainX ; i++) {
    for( j = 0 ; j < ScaleX ; j++) {
      MapX[icount++] = i;
    }
  }
  icount = 0;
  for( i = 0 ; i < domainY ; i++) {
    for( j = 0 ; j < ScaleY ; j++) {
      MapY[icount++] = i;
    }
  }
#endif // IPPL_GDL
}
//------------------------------------------------------------------
/***************************************************************************
 * $RCSfile: FieldView.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 * IPPL_VERSION_ID: $Id: FieldView.cpp,v 1.1.1.1 2003/01/23 07:40:33 adelmann Exp $ 
 ***************************************************************************/
