#include "LevelBld.H"

#include "Accel.H"
#include "Accel_F.H"
#include "Derive_F.H"

using std::string;

static Box the_same_box(const Box& b)
{
    return b;
}

static Box grow_box_by_one(const Box& b)
{
    return BoxLib::grow(b, 1);
}

typedef StateDescriptor::BndryFunc BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int norm_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_scalar_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i, scalar_bc[lo_bc[i]]);
        bc.setHi(i, scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, norm_vel_bc[lo_bc[0]]);
    bc.setHi(0, norm_vel_bc[hi_bc[0]]);
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    bc.setLo(1, norm_vel_bc[lo_bc[1]]);
    bc.setHi(1, norm_vel_bc[hi_bc[1]]);
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    bc.setLo(2, norm_vel_bc[lo_bc[2]]);
    bc.setHi(2, norm_vel_bc[hi_bc[2]]);
}

void
Accel::variable_setup()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();

    setup();

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    error_setup();
}

void
Accel::setup()
{

    // Initalize the parallel F90 constructs
    BL_FORT_PROC_CALL(PARALLEL_INIT_FORTRAN,parallel_init_fortran)();

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    BCRec bc;

    store_in_checkpoint = true;
    desc_lst.addDescriptor(Elec_Potential_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);
    store_in_checkpoint = false;
    desc_lst.addDescriptor(Elec_Field_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, BL_SPACEDIM,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);
    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(Elec_Potential_Type, 0, "Ephi", bc,
                          BndryFunc(BL_FORT_PROC_CALL(GENERIC_FILL,
                                                      generic_fill)));
    set_x_vel_bc(bc, phys_bc);
    desc_lst.setComponent(Elec_Field_Type, 0, "E_x", bc,
                          BndryFunc(BL_FORT_PROC_CALL(GENERIC_FILL,
                                                          generic_fill)));
   set_y_vel_bc(bc, phys_bc);
   desc_lst.setComponent(Elec_Field_Type, 1, "E_y", bc,
                         BndryFunc(BL_FORT_PROC_CALL(GENERIC_FILL,
                                                         generic_fill)));
   set_z_vel_bc(bc, phys_bc);
   desc_lst.setComponent(Elec_Field_Type, 2, "E_z", bc,
                         BndryFunc(BL_FORT_PROC_CALL(GENERIC_FILL,
                                                         generic_fill)));

   derive_lst.add("maggrav", IndexType::TheCellType(), 1,
                  BL_FORT_PROC_CALL(CA_DERMAGGRAV, ca_dermaggrav),
                  the_same_box);
   derive_lst.addComponent("maggrav", desc_lst, Elec_Field_Type, 0, BL_SPACEDIM);

    //
    // We want a derived type that corresponds to the number of particles in
    // each cell. We only intend to use it in plotfiles for debugging purposes.
    // We'll just use the DERNULL since don't do anything in fortran for now.
    // We'll actually set the values in `writePlotFile()`.
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, Elec_Field_Type, 0, 1);

    derive_lst.add("particle_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), grow_box_by_one);
    derive_lst.addComponent("particle_mass_density", desc_lst, Elec_Field_Type, 0, 1);

    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), the_same_box);
    derive_lst.addComponent("total_particle_count", desc_lst, Elec_Field_Type, 0, 1);

    derive_lst.add("total_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), grow_box_by_one);
    derive_lst.addComponent("total_density", desc_lst, Elec_Field_Type, 0, 1);
}

