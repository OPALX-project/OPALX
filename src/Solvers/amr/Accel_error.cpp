#include "Accel.H"
#include "Accel_error_F.H"

void
Accel::error_setup()
{
    err_list.add("total_particle_count", 1, ErrorRec::Standard,
                 BL_FORT_PROC_CALL(TAG_PART_CNT_ERR, tag_part_cnt_err));
}

void
Accel::manual_tags_placement (TagBoxArray&    tags,
                            Array<IntVect>& bf_lev)
{

}
