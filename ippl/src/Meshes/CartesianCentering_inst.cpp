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
#include "Meshes/CartesianCentering.h"

//-----------------------------------------------------------------------------
// Set up the CenteringEnum arrays describing these common cartesian centerings
// N.B.: the name "CCCEnums" is a shortened form of the original name for this
// class, "CommonCartesianCenteringEnums"
//
// Storage order of these CenteringEnum arrays:
// Where there is more than one component, the component index varies fastest,
// and the dimension index varies slowest--like a 2D Fortran-order array 
// dimensioned array(NComponents,D). With multicomponent types such as Tenzor,
// which intrinsically have 2 component indices (i,j), these are linearized 
// with the 1st component index i varying fastest and the 2nd component index j
// varying slowest. One-dimensional symmetric tensors (SymTenzor) have one 
// linearized component, 2D SymTenzor's have 3 linearized components, and 3D
// SymTenzor's have 6 linearized components; other components of SymTenzor's 
// are by definition zero and are not stored in the SymTenzor type, and so no 
// centering enum values are stored for them.
//-----------------------------------------------------------------------------

//11111111111111111111111111111111111111111111111111111111111111111111111111111
// 1D fields
//11111111111111111111111111111111111111111111111111111111111111111111111111111

// 1D field of scalars (or 1D vectors, or 1D tensors, or 1D sym. tensors)
// All components cell-centered:
CenteringEnum CCCEnums<1U,1U,0U>::allCell[1U*1U] =
{CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<1U,1U,0U>::allVertex[1U*1U] = 
{VERTEX};
// Componentwise centering along/perpendicular to component direction:
CenteringEnum CCCEnums<1U,1U,0U>::allFace[1U*1U] = 
{VERTEX};
CenteringEnum CCCEnums<1U,1U,0U>::allEdge[1U*1U] = 
{CELL};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<1U,1U,0U>::vectorFace[1U*1U] = 
{VERTEX};
CenteringEnum CCCEnums<1U,1U,0U>::vectorEdge[1U*1U] = 
{CELL};


//22222222222222222222222222222222222222222222222222222222222222222222222222222
// 2D fields
//22222222222222222222222222222222222222222222222222222222222222222222222222222

// 2D field of scalars (or 1D vectors, or 1D tensors, or 1D sym. tensors)
// All components cell-centered:
CenteringEnum CCCEnums<2U,1U,0U>::allCell[2U*1U] = 
{CELL, 
 CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<2U,1U,0U>::allVertex[2U*1U] = 
{VERTEX, 
 VERTEX};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<2U,1U,0U>::allFace[2U*1U] = 
{VERTEX, 
 CELL};
CenteringEnum CCCEnums<2U,1U,0U>::allEdge[2U*1U] = 
{CELL, 
 VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<2U,1U,1U>::allFace[2U*1U] = 
{CELL, 
 VERTEX};
CenteringEnum CCCEnums<2U,1U,1U>::allEdge[2U*1U] = 
{VERTEX, 
 CELL};

// 2D field of 2D vectors:
// All components cell-centered:
CenteringEnum CCCEnums<2U,2U,0U>::allCell[2U*2U] = 
{CELL, CELL, 
 CELL, CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<2U,2U,0U>::allVertex[2U*2U] = 
{VERTEX, VERTEX,
 VERTEX, VERTEX};
// Componentwise centering along/perpendicular to component direction:
CenteringEnum CCCEnums<2U,2U,0U>::vectorFace[2U*2U] = 
{VERTEX, CELL,
 CELL, VERTEX};
CenteringEnum CCCEnums<2U,2U,0U>::vectorEdge[2U*2U] = 
{CELL, VERTEX,
 VERTEX, CELL};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<2U,2U,0U>::allFace[2U*2U] = 
{VERTEX, VERTEX,
 CELL, CELL};
CenteringEnum CCCEnums<2U,2U,0U>::allEdge[2U*2U] = 
{CELL, CELL,
 VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<2U,2U,1U>::allFace[2U*2U] = 
{CELL, CELL,
 VERTEX, VERTEX};
CenteringEnum CCCEnums<2U,2U,1U>::allEdge[2U*2U] = 
{VERTEX, VERTEX,
 CELL, CELL};

// 2D field of 2D tensors:
// All components cell-centered:
CenteringEnum CCCEnums<2U,4U,0U>::allCell[2U*4U] = 
{CELL, CELL, CELL, CELL, 
 CELL, CELL, CELL, CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<2U,4U,0U>::allVertex[2U*4U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, 
 VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<2U,4U,0U>::allFace[2U*4U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, 
 CELL, CELL, CELL, CELL};
CenteringEnum CCCEnums<2U,4U,0U>::allEdge[2U*4U] = 
{CELL, CELL, CELL, CELL, 
 VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<2U,4U,1U>::allFace[2U*4U] = 
{CELL, CELL, CELL, CELL, 
 VERTEX, VERTEX, VERTEX, VERTEX};
CenteringEnum CCCEnums<2U,4U,1U>::allEdge[2U*4U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, 
 CELL, CELL, CELL, CELL};

// 2D field of 2D symmetric tensors:
// All components cell-centered:
CenteringEnum CCCEnums<2U,3U,0U>::allCell[2U*3U] = 
{CELL, CELL, CELL, 
 CELL, CELL, CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<2U,3U,0U>::allVertex[2U*3U] = 
{VERTEX, VERTEX, VERTEX, 
 VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<2U,3U,0U>::allFace[2U*3U] = 
{VERTEX, VERTEX, VERTEX, 
 CELL, CELL, CELL};
CenteringEnum CCCEnums<2U,3U,0U>::allEdge[2U*3U] = 
{CELL, CELL, CELL, 
 VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<2U,3U,1U>::allFace[2U*3U] = 
{CELL, CELL, CELL, 
 VERTEX, VERTEX, VERTEX};
CenteringEnum CCCEnums<2U,3U,1U>::allEdge[2U*3U] = 
{VERTEX, VERTEX, VERTEX, 
 CELL, CELL, CELL};


//33333333333333333333333333333333333333333333333333333333333333333333333333333
// 3D fields
//33333333333333333333333333333333333333333333333333333333333333333333333333333

// 3D field of scalars (or 1D vectors, or 1D tensors, or 1D sym. tensors)
// All components cell-centered:
CenteringEnum CCCEnums<3U,1U,0U>::allCell[3U*1U] = 
{CELL, 
 CELL, 
 CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<3U,1U,0U>::allVertex[3U*1U] = 
{VERTEX, 
 VERTEX, 
 VERTEX};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<3U,1U,0U>::allFace[3U*1U] = 
{VERTEX, 
 CELL, 
 CELL};
CenteringEnum CCCEnums<3U,1U,0U>::allEdge[3U*1U] = 
{CELL, 
 VERTEX, 
 VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<3U,1U,1U>::allFace[3U*1U] = 
{CELL, 
 VERTEX, 
 CELL};
CenteringEnum CCCEnums<3U,1U,1U>::allEdge[3U*1U] = 
{VERTEX, 
 CELL, 
 VERTEX};
// Face/Edge centering perpendicular to/along direction 2:
CenteringEnum CCCEnums<3U,1U,2U>::allFace[3U*1U] = 
{CELL, 
 CELL, 
 VERTEX};
CenteringEnum CCCEnums<3U,1U,2U>::allEdge[3U*1U] = 
{VERTEX, 
 VERTEX, 
 CELL};

// 3D field of 2D vectors:
// All components cell-centered:
CenteringEnum CCCEnums<3U,2U,0U>::allCell[3U*2U] = 
{CELL, CELL, 
 CELL, CELL, 
 CELL, CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<3U,2U,0U>::allVertex[3U*2U] = 
{VERTEX, VERTEX, 
 VERTEX, VERTEX, 
 VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<3U,2U,0U>::allFace[3U*2U] = 
{VERTEX, VERTEX, 
 CELL, CELL, 
 CELL, CELL};
CenteringEnum CCCEnums<3U,2U,0U>::allEdge[3U*2U] = 
{CELL, CELL, 
 VERTEX, VERTEX, 
 VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<3U,2U,1U>::allFace[3U*2U] = 
{CELL, CELL, 
 VERTEX, VERTEX, 
 CELL, CELL};
CenteringEnum CCCEnums<3U,2U,1U>::allEdge[3U*2U] = 
{VERTEX, VERTEX, 
 CELL, CELL, 
 VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 2:
CenteringEnum CCCEnums<3U,2U,2U>::allFace[3U*2U] = 
{CELL, CELL, 
 CELL, CELL, 
 VERTEX, VERTEX};
CenteringEnum CCCEnums<3U,2U,2U>::allEdge[3U*2U] = 
{VERTEX, VERTEX, 
 VERTEX, VERTEX, 
 CELL, CELL};

// 3D field of 3D vectors:
// All components cell-centered:
CenteringEnum CCCEnums<3U,3U,0U>::allCell[3U*3U] = 
{CELL, CELL, CELL, 
 CELL, CELL, CELL, 
 CELL, CELL, CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<3U,3U,0U>::allVertex[3U*3U] = 
{VERTEX, VERTEX, VERTEX, 
 VERTEX, VERTEX, VERTEX, 
 VERTEX, VERTEX, VERTEX};
// Componentwise centering along/perpendicular to component direction:
CenteringEnum CCCEnums<3U,3U,0U>::vectorFace[3U*3U] = 
{VERTEX, CELL, CELL,
 CELL, VERTEX, CELL,
 CELL, CELL, VERTEX};
CenteringEnum CCCEnums<3U,3U,0U>::vectorEdge[3U*3U] = 
{CELL, VERTEX, VERTEX,
 VERTEX, CELL, VERTEX,
 VERTEX, VERTEX, CELL};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<3U,3U,0U>::allFace[3U*3U] = 
{VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL,
 CELL, CELL, CELL};
CenteringEnum CCCEnums<3U,3U,0U>::allEdge[3U*3U] = 
{CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<3U,3U,1U>::allFace[3U*3U] = 
{CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL};
CenteringEnum CCCEnums<3U,3U,1U>::allEdge[3U*3U] = 
{VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 2:
CenteringEnum CCCEnums<3U,3U,2U>::allFace[3U*3U] = 
{CELL, CELL, CELL,
 CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX};
CenteringEnum CCCEnums<3U,3U,2U>::allEdge[3U*3U] = 
{VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL};

// 3D field of 3D tensors:
// All components cell-centered:
CenteringEnum CCCEnums<3U,9U,0U>::allCell[3U*9U] = 
{CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<3U,9U,0U>::allVertex[3U*9U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<3U,9U,0U>::allFace[3U*9U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL};
CenteringEnum CCCEnums<3U,9U,0U>::allEdge[3U*9U] = 
{CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<3U,9U,1U>::allFace[3U*9U] = 
{CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL};
CenteringEnum CCCEnums<3U,9U,1U>::allEdge[3U*9U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 2:
CenteringEnum CCCEnums<3U,9U,2U>::allFace[3U*9U] = 
{CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
CenteringEnum CCCEnums<3U,9U,2U>::allEdge[3U*9U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL, CELL};

// 3D field of 3D symmetric tensors:
// All components cell-centered:
CenteringEnum CCCEnums<3U,6U,0U>::allCell[3U*6U] = 
{CELL, CELL, CELL, CELL, CELL, CELL, 
 CELL, CELL, CELL, CELL, CELL, CELL, 
 CELL, CELL, CELL, CELL, CELL, CELL};
// All components vertex-centered:
CenteringEnum CCCEnums<3U,6U,0U>::allVertex[3U*6U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, 
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, 
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 0:
CenteringEnum CCCEnums<3U,6U,0U>::allFace[3U*6U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL,
 CELL, CELL, CELL, CELL, CELL, CELL};
CenteringEnum CCCEnums<3U,6U,0U>::allEdge[3U*6U] = 
{CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 1:
CenteringEnum CCCEnums<3U,6U,1U>::allFace[3U*6U] = 
{CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL};
CenteringEnum CCCEnums<3U,6U,1U>::allEdge[3U*6U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
// Face/Edge centering perpendicular to/along direction 2:
CenteringEnum CCCEnums<3U,6U,2U>::allFace[3U*6U] = 
{CELL, CELL, CELL, CELL, CELL, CELL,
 CELL, CELL, CELL, CELL, CELL, CELL,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX};
CenteringEnum CCCEnums<3U,6U,2U>::allEdge[3U*6U] = 
{VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 VERTEX, VERTEX, VERTEX, VERTEX, VERTEX, VERTEX,
 CELL, CELL, CELL, CELL, CELL, CELL};

// Names for the centering classes:
// allEdge scalar:
template<>   char* 
CartesianCentering<CCCEnums<2U,1U,0U>::allEdge,1U,1U>::
CenteringName = "CartesianCentering: allEdge centering";
template<>   char* 
CartesianCentering<CCCEnums<2U,1U,1U>::allEdge,1U,1U>::
CenteringName = "CartesianCentering: allEdge centering";
// allCell scalar:
template<>   char* 
CartesianCentering<CCCEnums<1U,1U,0U>::allCell,1U,0U>::
CenteringName = "CartesianCentering: allCell(1U,1U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<2U,2U,0U>::allCell,1U,0U>::
CenteringName = "CartesianCentering: allCell(2U,1U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<3U,3U,0U>::allCell,1U,0U>::
CenteringName = "CartesianCentering: allCell(3U,1U,0U) centering";
// allCell vector:
template<>   char* 
CartesianCentering<CCCEnums<2U,2U,0U>::allCell,2U,0U>::
CenteringName = "CartesianCentering: allCell(2U,2U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<3U,3U,0U>::allCell,3U,0U>::
CenteringName = "CartesianCentering: allCell(3U,3U,0U) centering";
// allVertex scalar:
template<>   char* 
CartesianCentering<CCCEnums<1U,1U,0U>::allVertex,1U,0U>::
CenteringName = "CartesianCentering: allVertex(1U,1U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<2U,2U,0U>::allVertex,1U,0U>::
CenteringName = "CartesianCentering: allVertex(2U,1U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<3U,3U,0U>::allVertex,1U,0U>::
CenteringName = "CartesianCentering: allVertex(3U,1U,0U) centering";
// allVertex vector:
template<>   char* 
CartesianCentering<CCCEnums<2U,2U,0U>::allVertex,2U,0U>::
CenteringName = "CartesianCentering: allVertex(2U,2U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<3U,3U,0U>::allVertex,3U,0U>::
CenteringName = "CartesianCentering: allVertex(3U,3U,0U) centering";
// vectorFace:
template<>   char* 
CartesianCentering<CCCEnums<1U,1U,0U>::vectorFace,1U,0U>::
CenteringName = "CartesianCentering: vectorFace(1U,1U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<2U,2U,0U>::vectorFace,2U,0U>::
CenteringName = "CartesianCentering: vectorFace(2U,2U,0U) centering";
template<>   char* 
CartesianCentering<CCCEnums<3U,3U,0U>::vectorFace,3U,0U>::
CenteringName = "CartesianCentering: vectorFace(3U,3U,0U) centering";
//.. fill these in more later (tjw) 

/***************************************************************************
 * $RCSfile: CartesianCentering_inst.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:28 $
 * IPPL_VERSION_ID: $Id: CartesianCentering_inst.cpp,v 1.1.1.1 2003/01/23 07:40:28 adelmann Exp $ 
 ***************************************************************************/
