// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef TAGS_H
#define TAGS_H

/*
 * Tags.h - list of special tags used by each major component in the ippl
 *	library.  When a new general communication cycle (i.e. swapping
 *	boundaries) is added, a new item should be added to this list.
 * BH, 6/19/95
 *
 * Updated for beta, JVWR, 7/27/96
 */

// special tag used to indicate the program should quit.  The values are
// arbitrary, but non-zero.
#define IPPL_ABORT_TAG            5    // program should abort()
#define IPPL_EXIT_TAG             6    // program should exit()
#define IPPL_RETRANSMIT_TAG       7    // node should resend a message
#define IPPL_MSG_OK_TAG           8    // some messages were sent correctly


// tags for reduction
#define COMM_REDUCE_SEND_TAG    1000000
#define COMM_REDUCE_RECV_TAG    1100000
#define COMM_REDUCE_SCATTER_TAG 1200000
#define COMM_REDUCE_CYCLE        100000


// tag for applying parallel periodic boundary condition.

#define BC_PARALLEL_PERIODIC_TAG 1500000
#define BC_TAG_CYCLE              100000

// Field<T,Dim> tags
#define F_GUARD_CELLS_TAG       2000000 // Field::fillGuardCells()
#define F_WRITE_TAG             2100000 // Field::write()
#define F_READ_TAG              2200000 // Field::read()
#define F_GEN_ASSIGN_TAG        2300000 // assign(BareField,BareField)
#define F_REPARTITION_BCAST_TAG 2400000 // broadcast in FieldLayout::repartion.
#define F_REDUCE_PERP_TAG       2500000 // reduction in binary load balance.
#define F_GETSINGLE_TAG         2600000 // IndexedBareField::getsingle()
#define F_REDUCE_TAG            2700000 // Reduction in minloc/maxloc
#define F_LAYOUT_IO_TAG         2800000 // Reduction in minloc/maxloc
#define F_TAG_CYCLE              100000

// Tags for FieldView and FieldBlock
#define FV_2D_TAG               3000000 // FieldView::update_2D_data()
#define FV_3D_TAG               3100000 // FieldView::update_2D_data()
#define FV_TAG_CYCLE             100000

#define FB_WRITE_TAG            3200000 // FieldBlock::write()
#define FB_READ_TAG             3300000 // FieldBlock::read()
#define FB_TAG_CYCLE             100000

#define FP_GATHER_TAG           3400000 // FieldPrint::print()
#define FP_TAG_CYCLE             100000

// Tags for DiskField
#define DF_MAKE_HOST_MAP_TAG    3500000
#define DF_FIND_RECV_NODES_TAG  3600000
#define DF_QUERY_TAG            3700000
#define DF_READ_TAG             3800000
#define DF_OFFSET_TAG           3900000
#define DF_READ_META_TAG        4000000
#define DF_TAG_CYCLE             100000

// Special tags used by Particle classes for communication.
#define P_WEIGHTED_LAYOUT_TAG   5000000
#define P_WEIGHTED_RETURN_TAG   5100000
#define P_WEIGHTED_TRANSFER_TAG 5200000
#define P_SPATIAL_LAYOUT_TAG    5300000
#define P_SPATIAL_RETURN_TAG    5400000
#define P_SPATIAL_TRANSFER_TAG  5500000
#define P_SPATIAL_GHOST_TAG     5600000
#define P_SPATIAL_RANGE_TAG     5700000
#define P_RESET_ID_TAG          5800000
#define P_LAYOUT_CYCLE           100000

// Tags for Ippl setup
#define IPPL_MAKE_HOST_MAP_TAG    6000000
#define IPPL_TAG_CYCLE             100000

// Tags for Conejo load balancer
#define F_CONEJO_BALANCER_TAG      7000000
#define F_CB_BCAST_TAG             7100000
#define F_CB_DOMAIN_TAG            7200000

// Tags for VnodeMultiBalancer
#define VNMB_PARTIAL_TAG           8000000
#define VNMB_COMPLETE_TAG          8100000
#define VNMB_TAG_CYCLE              100000

// Tags for Ippl application codes
#define IPPL_APP_TAG0    9000000
#define IPPL_APP_TAG1    9100000
#define IPPL_APP_TAG2    9200000
#define IPPL_APP_TAG3    9300000
#define IPPL_APP_TAG4    9400000
#define IPPL_APP_TAG5    9500000
#define IPPL_APP_TAG6    9600000
#define IPPL_APP_TAG7    9700000
#define IPPL_APP_TAG8    9800000
#define IPPL_APP_TAG9    9900000
#define IPPL_APP_CYCLE    100000

#endif // TAGS_H
