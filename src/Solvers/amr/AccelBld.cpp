
#include <LevelBld.H>
#include "Accel.H"

class AccelBld
    :
    public LevelBld
{
    virtual void variable_setup();
    virtual void variable_cleanup();

    // hack copies for BoxLib overriding
    virtual void variableSetUp();
    virtual void variableCleanUp();

    virtual AmrLevel *operator() ();
    virtual AmrLevel *operator() (Amr& papa, int lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba, Real time);
};

AccelBld Accel_bld;

LevelBld*
get_level_bld()
{
    return &Accel_bld;
}

void
AccelBld::variable_setup()
{
    Accel::variable_setup();
}

void
AccelBld::variable_cleanup()
{
    Accel::variable_cleanup();
}

AmrLevel*
AccelBld::operator() ()
{
    return new Accel;
}

AmrLevel*
AccelBld::operator() (Amr& papa, int lev, const Geometry& level_geom,
                    const BoxArray& ba, Real time)
{
    return new Accel(papa, lev, level_geom, ba, time);
}

// override hacks, copies of above
LevelBld*
getLevelBld()
{
    return &Accel_bld;
}

void AccelBld::variableSetUp()
{
    Accel::variable_setup();
}
void AccelBld::variableCleanUp()
{
    Accel::variable_cleanup();
}
