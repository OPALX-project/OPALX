
! ::: -----------------------------------------------------------

      subroutine generic_fill(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3, &
                              domlo,domhi,delta,xlo,time,bc)

      implicit none
      include 'bc_types.fi'
      integer var_l1,var_l2,var_l3,var_h1,var_h2,var_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision var(var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)

      call filcc(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,domlo,domhi,delta,xlo,bc)

      end subroutine generic_fill

