
! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine PROBINIT (init,name,namlen,problo,probhi)

      implicit none

      integer init, namlen
      integer name(namlen)
      double precision problo(3), probhi(3)

      end

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine parallel_init_fortran()

        ! Passing data from C++ into f90
        use parallel

!         call parallel_initialize()

      end subroutine parallel_init_fortran
