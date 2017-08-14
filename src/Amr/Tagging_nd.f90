!> @file Tagging_nd.f90
!! @author Weiqun Zhang
!! @date October 2016, LBNL
!! @details This Fortran routine is used for tagging
!!          the cells at a level for refinement.
!!          In case of refinement a cell of the tagbox is
!!          marked "set" otherwise "clear".
!! @brief Fortran routine for tagging cells.

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the state
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag        <=  integer tag array
! ::: tag_lo,hi   => index extent of tag array
! ::: state       => state array
! ::: state_lo,hi => index extent of state array
! ::: set         => integer value to tag cell for refinement
! ::: clear       => integer value to untag cell
! ::: lo,hi       => work region we are allowed to change
! ::: dx          => cell size
! ::: problo      => phys loc of lower left corner of prob domain
! ::: time        => problem evolution time
! ::: level       => refinement level of this array
! ::: -----------------------------------------------------------

!> Called by AmrOpal::ErrorEst.
!! @param tag integer tag array
!! @param tag_lo lower index extent of tag array
!! @param tag_hi upper index extent of tag array
!! @param state is the state array
!! @param set is the integer value to tag cell fo refinement
!! @param clear is the integer value to untag a cell
!! @param lo is the lower left corner of the work region we are allowed to change
!! @param hi is the upper right corner of the work region we are allowed to change
!! @param dx is the cell size
!! @param problo is the physical location of the lower left corner of the problem domain
!! @param time is the problem evolution time
!! @param level is the refinement level of this array
subroutine state_error(tag,tag_lo,tag_hi, &
                       state,state_lo,state_hi, &
                       set,clear,&
                       lo,hi,&
                       dx,problo,time,phierr) bind(C, name="state_error")

    implicit none
    
    integer          :: lo(3),hi(3)
    integer          :: state_lo(3),state_hi(3)
    integer          :: tag_lo(3),tag_hi(3)
    double precision :: state(state_lo(1):state_hi(1), &
                              state_lo(2):state_hi(2), &
                              state_lo(3):state_hi(3))
    integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    double precision :: problo(3),dx(3),time,phierr
    integer          :: set,clear
    
    integer          :: i, j, k
    
    do       k = lo(3), hi(3)
        do    j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (abs(state(i,j,k)) .ge. phierr) then
                   tag(i,j,k) = set
                endif
            enddo
        enddo
    enddo
end subroutine state_error



subroutine tag_potential_strength(tag, tag_lo, tag_hi, &
                                  state, state_lo, state_hi, &
                                  set, clear, &
                                  lo, hi, &
                                  dx, problo, time, phi) bind(C, name="tag_potential_strength")
    
    implicit none
    
    integer          :: lo(3),hi(3)
    integer          :: state_lo(3),state_hi(3)
    integer          :: tag_lo(3),tag_hi(3)
    double precision :: state(state_lo(1):state_hi(1), &
                              state_lo(2):state_hi(2), &
                              state_lo(3):state_hi(3))
    integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    double precision :: problo(3),dx(3),time,phi
    integer          :: set,clear
    
    integer          :: i, j, k
    
    ! Tag on regions of high phi
    do       k = lo(3), hi(3)
        do    j = lo(2), hi(2)
            do i = lo(1), hi(1)
                if (abs(state(i,j,k)) .ge. phi) then
                    tag(i,j,k) = set
                endif
            enddo
        enddo
    enddo
    
end subroutine tag_potential_strength
