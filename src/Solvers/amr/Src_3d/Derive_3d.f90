
      subroutine ca_derlapvar(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,nv, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! Compute the weighted-Laplacian of the variable for use in error estimation
      !
      implicit none

      integer          lo(3), hi(3)
      integer          var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision var(var_l1:var_h1,var_l2:var_h2,var_l3:var_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      double precision ::  delu(3,var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)
      double precision :: delua(3,var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)
      double precision :: delu2(9), delu3(9), delu4(9)
      double precision :: num, denom
      integer          :: i,j,k

      ! This value is taken from FLASH
      double precision, parameter:: epsil=0.02d0

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      !$OMP PARALLEL PRIVATE(i,j,k)
      !$OMP DO
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(1,i,j,k) =     dat(i+1,j,k,1) -      dat(i-1,j,k,1) 
         delua(1,i,j,k) = abs(dat(i+1,j,k,1)) + abs(dat(i-1,j,k,1))
      end do
      end do
      end do
      !$OMP END DO NOWAIT

      ! d/dy
      !$OMP DO
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(2,i,j,k) =     dat(i,j+1,k,1)  -     dat(i,j-1,k,1) 
         delua(2,i,j,k) = abs(dat(i,j+1,k,1)) + abs(dat(i,j-1,k,1))
      end do
      end do
      end do
      !$OMP END DO NOWAIT

      ! d/dz
      !$OMP DO
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(3,i,j,k) =     dat(i,j,k+1,1)  -     dat(i,j,k-1,1)
         delua(3,i,j,k) = abs(dat(i,j,k+1,1)) + abs(dat(i,j,k-1,1))
      end do
      end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      !$OMP PARALLEL DO PRIVATE(i,j,k,delu2,delu3,delu4,num,denom)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         

         ! d/dxdx
          delu2(1) =     delu(1,i+1,j,k)  -     delu(1,i-1,j,k)
          delu3(1) = abs(delu(1,i+1,j,k)) + abs(delu(1,i-1,j,k))
          delu4(1) =    delua(1,i+1,j,k)  +    delua(1,i-1,j,k)

          ! d/dydx
          delu2(2) =     delu(1,i,j+1,k)  -     delu(1,i,j-1,k)
          delu3(2) = abs(delu(1,i,j+1,k)) + abs(delu(1,i,j-1,k))
          delu4(2) =    delua(1,i,j+1,k)  +    delua(1,i,j-1,k)

          ! d/dxdy
          delu2(3) =     delu(2,i+1,j,k)  -     delu(2,i-1,j,k)
          delu3(3) = abs(delu(2,i+1,j,k)) + abs(delu(2,i-1,j,k))
          delu4(3) =    delua(2,i+1,j,k)  +    delua(2,i-1,j,k)

          ! d/dydy
          delu2(4) =     delu(2,i,j+1,k)  -     delu(2,i,j-1,k)
          delu3(4) = abs(delu(2,i,j+1,k)) + abs(delu(2,i,j-1,k))
          delu4(4) =    delua(2,i,j+1,k)  +    delua(2,i,j-1,k)

          ! d/dzdx
          delu2(5) =     delu(1,i,j,k+1)  -     delu(1,i,j,k-1)
          delu3(5) = abs(delu(1,i,j,k+1)) + abs(delu(1,i,j,k-1))
          delu4(5) =    delua(1,i,j,k+1)  +    delua(1,i,j,k-1)

          ! d/dzdy
          delu2(6) =     delu(2,i,j,k+1)  -     delu(2,i,j,k-1)
          delu3(6) = abs(delu(2,i,j,k+1)) + abs(delu(2,i,j,k-1))
          delu4(6) =    delua(2,i,j,k+1)  +    delua(2,i,j,k-1)

          ! d/dxdz
          delu2(7) =     delu(3,i+1,j,k)  -     delu(3,i-1,j,k)
          delu3(7) = abs(delu(3,i+1,j,k)) + abs(delu(3,i-1,j,k))
          delu4(7) =    delua(3,i+1,j,k)  +    delua(3,i-1,j,k)

          ! d/dydz
          delu2(8) =     delu(3,i,j+1,k)  -     delu(3,i,j-1,k)
          delu3(8) = abs(delu(3,i,j+1,k)) + abs(delu(3,i,j-1,k))
          delu4(8) =    delua(3,i,j+1,k)  +    delua(3,i,j-1,k)

          ! d/dzdz
          delu2(9) =     delu(3,i,j,k+1)  -     delu(3,i,j,k-1)
          delu3(9) = abs(delu(3,i,j,k+1)) + abs(delu(3,i,j,k-1))
          delu4(9) =    delua(3,i,j,k+1)  +    delua(3,i,j,k-1)

         ! compute the error
         num   =  delu2(1)**2 + delu2(2)**2 + delu2(3)**2 + delu2(4)**2 &
                 +delu2(5)**2 + delu2(6)**2 + delu2(7)**2 + delu2(8)**2 &
                 +delu2(9)**2

         denom = (delu3(1) + (epsil*delu4(1)+1.d-99))**2 + &
                 (delu3(2) + (epsil*delu4(2)+1.d-99))**2 + &
                 (delu3(3) + (epsil*delu4(3)+1.d-99))**2 + &
                 (delu3(4) + (epsil*delu4(4)+1.d-99))**2 + &
                 (delu3(5) + (epsil*delu4(5)+1.d-99))**2 + &
                 (delu3(6) + (epsil*delu4(6)+1.d-99))**2 + &
                 (delu3(7) + (epsil*delu4(7)+1.d-99))**2 + &
                 (delu3(8) + (epsil*delu4(8)+1.d-99))**2 + &
                 (delu3(9) + (epsil*delu4(9)+1.d-99))**2

         var(i,j,k,1) = sqrt(num/denom)

      end do
      end do
      end do
      !$OMP END PARALLEL DO


      end subroutine ca_derlapvar

!-----------------------------------------------------------------------

      subroutine ca_derstate(state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv, &
                             dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! The incoming   "dat" vector contains (rho,T,(rho X)_1)
      ! The outgoing "state" vector contains (rho,T,X_1)
      !
      implicit none 

      integer          lo(3), hi(3)
      integer          state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k

      if (nv .ne. 3) then
          print *,'... confusion in derstate ... nv should be 3 but is ',nv
          call bl_error('Error:: Derive_3d.f90 :: ca_derstate')
      end if
      !
      ! Density
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,1) = dat(i,j,k,1)
            end do
         end do
      end do
      !
      ! Temperature
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,2) = dat(i,j,k,2)
            end do
         end do
      end do
      !
      ! (rho X)_1 --> X_1
      !
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,3) = dat(i,j,k,3) / dat(i,j,k,1)
            end do
         end do
      end do
      !$OMP END PARALLEL DO
 
      end subroutine ca_derstate

!-----------------------------------------------------------------------

      subroutine ca_dermaggrav(maggrav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,ng, &
                               dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! Derive magnitude of the gravity vector.
      !
      implicit none 

      integer          lo(3), hi(3)
      integer          grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,ng
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision maggrav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3,ng)
      double precision     dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      ! 
      ! Here dat contains (grav_x, grav_y, grav_z)
      ! 
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               maggrav(i,j,k,1) = sqrt( dat(i,j,k,1)**2  + &
                                        dat(i,j,k,2)**2  + &
                                        dat(i,j,k,3)**2 )
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      end subroutine ca_dermaggrav

!-----------------------------------------------------------------------

      subroutine ca_dernull(kineng,ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine is used by particle_count.  Yes it does nothing.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision kineng(ken_l1:ken_h1,ken_l2:ken_h2,ken_l3:ken_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      end subroutine ca_dernull

!-----------------------------------------------------------------------
