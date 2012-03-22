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
#include "Field/Field.h"
#include "Field/IndexedField.h"
#include "DataSource/MakeDataSource.h"
#include "Index/SIndex.h"
#include "SubField/SubField.h"
#include "Utility/IpplStats.h"
#include "Profile/Profiler.h"

//=============================================================================
// Global functions
//=============================================================================

//-----------------------------------------------------------------------------
// Wrappers for internal mesh object constructor calls:
//-----------------------------------------------------------------------------
// Generic case, report unimplemented mesh type:
// (Compiler limitations doesn't allow this yet; user will get an obscure
// compile error instead of this runtime error if he specifies an unsupported
// mesh type for the Field Mesh parameter --tjw):
// template<class T, unsigned Dim, class M, class C>
// M* makeMesh(Field<T,Dim,M,C>& f)
// {
//   ERRORMSG("makeMesh() invoked from Field(): unimplemented mesh type" << endl);
// }

// Generic makeMesh function
template<class T, unsigned Dim, class M, class C>
inline M* makeMesh(Field<T,Dim,M,C>& f)
{ 
  TAU_TYPE_STRING(taustr, "MeshType* (" + CT(f) + " )" );
  TAU_PROFILE("makeMesh()", taustr, TAU_FIELD);
  
  NDIndex<Dim> ndi;
  ndi = f.getLayout().getDomain();
  return (new M(ndi));
}

/*
// Specialization for UniformCartesian
template<class T, unsigned Dim, class MFLOAT, class C>
UniformCartesian<Dim,MFLOAT>*
makeMesh(Field<T,Dim,UniformCartesian<Dim,MFLOAT>,C>& f)
{ 
  TAU_TYPE_STRING(taustr, "UniformCartesian<Dim,MFLOAT>* (" + CT(f) + " )" );
  TAU_PROFILE("makeMesh()", taustr, TAU_FIELD);
  
  NDIndex<Dim> ndi;
  ndi = f.getLayout().getDomain();
  return (new UniformCartesian<Dim,MFLOAT>(ndi));
}

// Specialization for Cartesian
template<class T, unsigned Dim, class MFLOAT, class C>
Cartesian<Dim,MFLOAT>*
makeMesh(Field<T,Dim,Cartesian<Dim,MFLOAT>,C>& f)
{
  TAU_TYPE_STRING(taustr, "Cartesian<Dim,MFLOAT>* (" + CT(f) + " )" );
  TAU_PROFILE("makeMesh()", taustr, TAU_FIELD);

  NDIndex<Dim> ndi;
  ndi = f.getLayout().getDomain();
  return (new Cartesian<Dim,MFLOAT>(ndi));
}
*/

//=============================================================================
// Field member functions
//=============================================================================


//////////////////////////////////////////////////////////////////////////
// A default constructor, which should be used only if the user calls the
// 'initialize' function before doing anything else.  There are no special
// checks in the rest of the Field methods to check that the Field has
// been properly initialized
template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field() {
  TAU_TYPE_STRING(taustr, CT(*this) + " void ()" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(0, true);
}


//////////////////////////////////////////////////////////////////////////
// Field destructor
template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::~Field() {
  TAU_TYPE_STRING(taustr, CT(*this) + " void ()" );
  TAU_PROFILE("Field::~Field()", taustr, TAU_FIELD); 
  delete_mesh();
}


//////////////////////////////////////////////////////////////////////////
// Create a new Field with a given layout and optional guard cells.
// The default type of BCond lets you add new ones dynamically.
// The makeMesh() global function is a way to allow for different types of
// constructor arguments for different mesh types.
template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Layout_t & l) : BareField<T,Dim>(l) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (" + CT(l) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Layout_t & l, const GuardCellSizes<Dim>& gc)
  : BareField<T,Dim>(l,gc) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (" + CT(l) + ", " + CT(gc) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Layout_t & l, const BConds<T,Dim,M,C>& bc)
  : BareField<T,Dim>(l), Bc(bc) {
  TAU_TYPE_STRING(taustr, "void (" + CT(l) + ", " + CT(bc) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Layout_t & l, const GuardCellSizes<Dim>& gc,
			const BConds<T,Dim,M,C>& bc)
  : BareField<T,Dim>(l,gc), Bc(bc) {
  TAU_TYPE_STRING(taustr, "void (" + CT(l) + ", " + CT(gc) + ", " + CT(bc) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Layout_t & l, const BConds<T,Dim,M,C>& bc,
			const GuardCellSizes<Dim>& gc)
  : BareField<T,Dim>(l,gc), Bc(bc) {
  TAU_TYPE_STRING(taustr, "void (" + CT(l) + ", " + CT(bc) + ", " + CT(gc) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(FieldSpec<T,Dim,M,C>& spec)
  : BareField<T,Dim>( (Layout_t &) spec.get_Layout(), spec.get_GC()),
    Bc(spec.get_BC()) {
  TAU_TYPE_STRING(taustr, "void (" + CT(spec) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(makeMesh(*this), true);
}


//////////////////////////////////////////////////////////////////////////
// Constructors which include a Mesh object as argument
template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Mesh_t& m, Layout_t & l)
  : BareField<T,Dim>(l) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (Mesh_t, " + CT(l) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Mesh_t& m, Layout_t & l,
			const GuardCellSizes<Dim>& gc)
  : BareField<T,Dim>(l,gc) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (Mesh_t, " + CT(l) + ", " + CT(gc) 
    + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Mesh_t& m, Layout_t & l,
			const BConds<T,Dim,M,C>& bc)
  : BareField<T,Dim>(l), Bc(bc) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(l) + ", " + CT(bc) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Mesh_t& m, Layout_t & l,
			const GuardCellSizes<Dim>& gc,
			const BConds<T,Dim,M,C>& bc)
  : BareField<T,Dim>(l,gc), Bc(bc) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(l) + ", " + CT(gc) + ", " +  CT(bc) 
    + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Mesh_t& m, Layout_t & l,
			const BConds<T,Dim,M,C>& bc,
			const GuardCellSizes<Dim>& gc)
  : BareField<T,Dim>(l,gc), Bc(bc) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(l) + ", " + CT(bc) + ", " +  CT(gc) 
    + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
Field<T,Dim,M,C>::Field(Mesh_t& m, FieldSpec<T,Dim,M,C>& spec)
  : BareField<T,Dim>( (Layout_t &) spec.get_Layout(), spec.get_GC()),
    Bc(spec.get_BC()) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(spec) + " )" );
  TAU_PROFILE("Field::Field()", taustr, TAU_FIELD); 
  store_mesh(&m, false);
}


//////////////////////////////////////////////////////////////////////////
// Initialize the Field, if it was constructed from the default constructor.
// This should NOT be called if the Field was constructed by providing
// a FieldLayout or FieldSpec
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Layout_t & l) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (" + CT(l) + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l);
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Layout_t & l,
				  const GuardCellSizes<Dim>& gc) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (" + CT(l) + ", " + CT(gc) +  " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l,gc);
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Layout_t & l,
				  const BConds<T,Dim,M,C>& bc) {
  TAU_TYPE_STRING(taustr, "void (" + CT(l) + ", " + CT(bc) +  " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l);
  Bc = bc;
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Layout_t & l,
				  const GuardCellSizes<Dim>& gc,
				  const BConds<T,Dim,M,C>& bc) {
  TAU_TYPE_STRING(taustr, "void (" + CT(l) + ", " + CT(gc) + ", " + CT(bc) +  " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l,gc);
  Bc = bc;
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Layout_t & l,
				  const BConds<T,Dim,M,C>& bc,
				  const GuardCellSizes<Dim>& gc) {
  TAU_TYPE_STRING(taustr, "void (" + CT(l) + ", " + CT(bc) + ", " + CT(gc) +  " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l,gc);
  Bc = bc;
  store_mesh(makeMesh(*this), true);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(FieldSpec<T,Dim,M,C>& spec) {
  TAU_TYPE_STRING(taustr, "void (" + CT(spec) + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize( (Layout_t &) spec.get_Layout(),
				spec.get_GC());
  Bc = spec.get_BC();
  store_mesh(makeMesh(*this), true);
}


//////////////////////////////////////////////////////////////////////////
// Initialize the Field, also specifying a mesh
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Mesh_t& m, Layout_t & l) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (Mesh_t, " + CT(l) + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l);
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Mesh_t& m, Layout_t & l,
				  const GuardCellSizes<Dim>& gc) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (Mesh_t, " + CT(l) + ", " + CT(gc) 
    + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l,gc);
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Mesh_t& m, Layout_t & l,
				  const BConds<T,Dim,M,C>& bc) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(l) + ", " + CT(bc) + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l);
  Bc = bc;
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Mesh_t& m, Layout_t & l,
				  const GuardCellSizes<Dim>& gc,
				  const BConds<T,Dim,M,C>& bc) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(l) + ", " + CT(gc) + ", " +  CT(bc) 
    + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l,gc);
  Bc = bc;
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Mesh_t& m, Layout_t & l,
				  const BConds<T,Dim,M,C>& bc,
				  const GuardCellSizes<Dim>& gc) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(l) + ", " + CT(bc) + ", " +  CT(gc) 
    + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize(l,gc);
  Bc = bc;
  store_mesh(&m, false);
}

template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::initialize(Mesh_t& m, FieldSpec<T,Dim,M,C>& spec) {
  TAU_TYPE_STRING(taustr, "void (Mesh_t, " + CT(spec) + " )" );
  TAU_PROFILE("Field::initialize()", taustr, TAU_FIELD); 
  BareField<T,Dim>::initialize( (Layout_t &) spec.get_Layout(),
				spec.get_GC());
  Bc = spec.get_BC();
  store_mesh(&m, false);
}


//////////////////////////////////////////////////////////////////////////
// If you make any modifications using an iterator, you must call this.
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::fillGuardCells(bool reallyFill) const {
  TAU_TYPE_STRING(taustr, CT(*this) + " void ()" );
  TAU_PROFILE("Field::fillGuardCells()", taustr, TAU_FIELD); 

  // Fill the internal guard cells.
  BareField<T,Dim>::fillGuardCells(reallyFill);

  // Handle the user-supplied boundary conditions. If we're not supposed
  // to really fill the guard cells and the BC does not change physical
  // cells, don't bother.
  if (reallyFill || Bc.changesPhysicalCells()) {
    // cast away const, so we can apply BC's
    Field<T,Dim,M,C>& ncf = const_cast<Field<T,Dim,M,C>&>(*this);
    ncf.getBConds().apply(ncf);
    INCIPPLSTAT(incBoundaryConditions);
  }
}


//////////////////////////////////////////////////////////////////////////
// When we apply a bracket it converts the type
// to IndexedField so that we can check at compile time
// that we have the right number of indexes and brackets.
template<class T, unsigned Dim, class M, class C>
IndexedField<T,Dim,1,M,C> Field<T,Dim,M,C>::operator[](const Index& idx) {
  TAU_TYPE_STRING(taustr, CT(*this) + " IndexedField<T,Dim,1,M,C> (Index )" );
  TAU_PROFILE("Field::operator[]()", taustr, TAU_FIELD); 
  return IndexedField<T,Dim,1,M,C>(*this,idx);
}


//////////////////////////////////////////////////////////////////////////
// Also allow using an integer instead of a whole Index
template<class T, unsigned Dim, class M, class C>
IndexedField<T,Dim,1,M,C> Field<T,Dim,M,C>::operator[](int i) {
  TAU_TYPE_STRING(taustr, CT(*this) + " IndexedField<T,Dim,1,M,C> (int )" );
  TAU_PROFILE("Field::operator[]()", taustr, TAU_FIELD); 
  return IndexedField<T,Dim,1,M,C>(*this,i);
}


//////////////////////////////////////////////////////////////////////////
// Also allow using a single NDIndex instead of N Index objects:
template<class T, unsigned D, class M, class C>
IndexedField<T,D,D,M,C> Field<T,D,M,C>::operator[](const NDIndex<D>& n) {
  TAU_TYPE_STRING(taustr, CT(*this) + " IndexedField<T,D,D,M,C> (" + CT(n) + " )" );
  TAU_PROFILE("Field::operator[]()", taustr, TAU_FIELD); 
  return IndexedField<T,D,D,M,C>(*this,n);
}


//////////////////////////////////////////////////////////////////////////
// Also allow using a sparse index
template<class T, unsigned D, class M, class C>
SubField<T,D,M,C,SIndex<D> > Field<T,D,M,C>::operator[](const SIndex<D>& s) {
  TAU_TYPE_STRING(taustr, CT(*this) + " SubField<T,D,M,C,SIndex<D>> (" +CT(s) + " )" );
  TAU_PROFILE("Field::operator[]()", taustr, TAU_FIELD); 
  return SubField<T,D,M,C,SIndex<D> >(*this, s);
}


//////////////////////////////////////////////////////////////////////////
// I/O (special stuff not inherited from BareField):
// Print out contents of Centering class
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::print_Centerings(std::ostream& out) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (ostream )" );
  TAU_PROFILE("Field::print_Centerings()", taustr, TAU_FIELD | TAU_IO); 
  Centering_t::print_Centerings(out);
}


//////////////////////////////////////////////////////////////////////////
// Repartition onto a new layout, or when the mesh changes
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::Repartition(UserList *userlist) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (UserList * )" );
  TAU_PROFILE("Field::Repartition()", taustr, TAU_FIELD ); 

  // see if the userlist corresponds to our current mesh
  if (mesh != 0 && mesh->get_Id() == userlist->getUserListID()) {
    // for now, we don't care if the mesh changes ... but we might later
  } else {
    // the layout has changed, so redistribute our data

    // Cast to the proper type of FieldLayout.
    Layout_t *newLayout = (Layout_t *)( userlist );

    // Build a new temporary field on the new layout.
    Field<T,Dim,M,C> tempField(get_mesh(), *newLayout, this->getGC(), getBConds());

    // Copy our data over to the new layout.
    tempField = *this;

    // Copy back the pointers to the new local fields.
    BareField<T,Dim>::Locals_ac = tempField.Locals_ac;

    INCIPPLSTAT(incRepartitions);
  }
}


//////////////////////////////////////////////////////////////////////////
// Tell this object that an object is being deleted
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::notifyUserOfDelete(UserList *userlist) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (UserList * )" );
  TAU_PROFILE("Field::notifyUserOfDelete()", taustr, TAU_FIELD ); 
  // see if the userlist corresponds to our current mesh
  if (mesh != 0 && mesh->get_Id() == userlist->getUserListID()) {
    mesh = 0;
  } else {
    // since this is not for our mesh, defer to the base class function
    BareField<T,Dim>::notifyUserOfDelete(userlist);
  }
}


//////////////////////////////////////////////////////////////////////////
// a virtual function which is called by this base class to get a
// specific instance of DataSourceObject based on the type of data
// and the connection method (the argument to the call).
template<class T, unsigned Dim, class M, class C>
DataSourceObject *Field<T,Dim,M,C>::createDataSourceObject(const char *nm,
							   DataConnect *dc,
							   int tm) {
  TAU_TYPE_STRING(taustr, CT(*this) 
    + "DataSourceObject * (char *, DataConnect *, int)");
  TAU_PROFILE("Field::createDataSourceObject()", taustr, TAU_FIELD ); 
  return make_DataSourceObject(nm, dc, tm, *this);
}


//////////////////////////////////////////////////////////////////////////
// store the given mesh object pointer, and the flag whether we use it or not
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::store_mesh(Mesh_t* m, bool WeOwn) {
  TAU_TYPE_STRING(taustr, CT(*this) + " void (Mesh_t *, bool)");
  TAU_PROFILE("Field::store_mesh()", taustr, TAU_FIELD ); 
  mesh = m;
  WeOwnMesh = WeOwn;
  if (mesh != 0)
    mesh->checkin(*this);
}


//////////////////////////////////////////////////////////////////////////
// delete the mesh object, if necessary; otherwise, just zero the pointer
template<class T, unsigned Dim, class M, class C>
void Field<T,Dim,M,C>::delete_mesh() {
  TAU_TYPE_STRING(taustr, CT(*this) + " void ()");
  TAU_PROFILE("Field::delete_mesh()", taustr, TAU_FIELD ); 
  if (mesh != 0) {
    mesh->checkout(*this);
    if (WeOwnMesh)
      delete mesh;
    mesh = 0;
  }
}


/***************************************************************************
 * $RCSfile: Field.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: Field.cpp,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/
