! :: ----------------------------------------------------------
! :: Average the fine grid phi onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! :: Note this differs from fort_avgdown in that there is no volume weighting.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine fort_avgdown_phi (crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
                                   fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                                   lo,hi,lrat)
      implicit none
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer lo(3), hi(3)
      integer lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)

      integer i, j, k, ic, jc, kc, ioff, joff, koff
      integer lratx, lraty, lratz
      double precision volfrac

      lratx   = lrat(1)
      lraty   = lrat(2)
      lratz   = lrat(3)
      volfrac = 1.d0/float(lrat(1)*lrat(2)*lrat(3))
      !
      ! ::::: set coarse grid to zero on overlap
      !
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = 0.d0
            enddo
         enddo
      enddo
      !
      ! ::::: sum fine data
      !
      do koff = 0, lratz-1
        !$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,ioff,joff)
        do kc = lo(3), hi(3)
          k = kc*lratz + koff
          do joff = 0, lraty-1
            do jc = lo(2), hi(2)
              j = jc*lraty + joff
              do ioff = 0, lratx-1
                do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  crse(ic,jc,kc) = crse(ic,jc,kc) + fine(i,j,k)
                enddo
              enddo
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo

      !$OMP PARALLEL DO PRIVATE(ic,jc,kc)
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = volfrac*crse(ic,jc,kc)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine fort_avgdown_phi

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine fort_edge_interp(flo, fhi, nc, ratio, dir, &
           fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)

      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      double precision &
           fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,n,P,M,L
      double precision val, df
      !
      ! Do linear in dir, pc transverse to dir, leave alone the fine values
      ! lining up with coarse edges--assume these have been set to hold the 
      ! values you want to interpolate to the rest.
      !
      if (dir.eq.0) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                     df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(i+M,P,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1)-ratio(dir),ratio(1)
                  do i=flo(0),fhi(0)
                     df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(P,j+M,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=flo(2),fhi(2)-ratio(dir),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0),ratio(0)
                     df = fine(i,j,k+ratio(dir),n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                              fine(P,L,k+M,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      endif

      end subroutine fort_edge_interp

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine fort_pc_edge_interp(lo, hi, nc, ratio, dir, &
           crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2, &
           fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer lo(3),hi(3), nc, ratio(0:3-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      double precision &
           crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      double precision &
           fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P
      !
      ! For edge-based data, fill fine values with piecewise-constant interp of coarse data.
      ! Operate only on faces that overlap--ie, only fill the fine faces that make up each
      ! coarse face, leave the in-between faces alone.
      !
      if (dir.eq.0) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(1)-1
                           fine(ii,jj+L,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(1)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj+P,kk,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      end subroutine fort_pc_edge_interp

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine fort_avg_ec_to_cc(lo, hi, bc_lo, &
           symmetry_type, &
           cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
           ecx, ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3, &
           ecy, ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3, &
           ecz, eczl1, eczl2, eczl3, eczh1, eczh2, eczh3, &
           problo)

      implicit none
      integer          :: lo(3),hi(3)
      integer          :: bc_lo(3)
      integer          :: symmetry_type
      integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
      integer          :: ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3
      integer          :: ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3
      integer          :: eczl1, eczl2, eczl3, eczh1, eczh2, eczh3
      double precision :: cc(ccl1:cch1,ccl2:cch2,ccl3:cch3, 3)
      double precision :: ecx(ecxl1:ecxh1,ecxl2:ecxh2, ecxl3:ecxh3)
      double precision :: ecy(ecyl1:ecyh1,ecyl2:ecyh2, ecyl3:ecyh3)
      double precision :: ecz(eczl1:eczh1,eczl2:eczh2, eczl3:eczh3)
      double precision :: problo(3)

      ! Local variables
      integer          :: i,j,k

      !$OMP PARALLEL PRIVATE(i,j,k)
      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,1) = 0.5d0 * ( ecx(i+1,j,k) + ecx(i,j,k) )
            enddo
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,2) = 0.5d0 * ( ecy(i,j+1,k) + ecy(i,j,k) )
            enddo
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,3) = 0.5d0 * ( ecz(i,j,k+1) + ecz(i,j,k) )
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      ! Note this assumes the lo end of the domain is 0.
      !
      if (ccl1 .lt. 0 .and. bc_lo(1) .eq. symmetry_type) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),-1
                  cc(i,j,k,1) = -cc(-i-1,j,k,1)
                  cc(i,j,k,2) =  cc(-i-1,j,k,2)
                  cc(i,j,k,3) =  cc(-i-1,j,k,3)
               enddo
            enddo
         enddo
      endif
      !
      ! Note this assumes the lo end of the domain is 0.
      !
      if (ccl2 .lt. 0 .and. bc_lo(2) .eq. symmetry_type) then
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               do j=lo(2),-1
                  cc(i,j,k,1) =  cc(i,-j-1,k,1)
                  cc(i,j,k,2) = -cc(i,-j-1,k,2)
                  cc(i,j,k,3) =  cc(i,-j-1,k,3)
               enddo
            enddo
         enddo
      endif
      !
      ! Note this assumes the lo end of the domain is 0.
      !
      if (ccl3 .lt. 0 .and. bc_lo(3) .eq. symmetry_type) then
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do k=lo(3),-1
                  cc(i,j,k,1) =  cc(i,j,-k-1,1)
                  cc(i,j,k,2) =  cc(i,j,-k-1,2)
                  cc(i,j,k,3) = -cc(i,j,-k-1,3)
               enddo
            enddo
         enddo
      endif
         
      end subroutine fort_avg_ec_to_cc

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine fort_compute_ec(lo, hi, &
           cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
           ecx, ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3, &
           ecy, ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3, &
           ecz, eczl1, eczl2, eczl3, eczh1, eczh2, eczh3, &
           dx)

      implicit none
      integer          :: lo(3),hi(3)
      integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
      integer          :: ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3
      integer          :: ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3
      integer          :: eczl1, eczl2, eczl3, eczh1, eczh2, eczh3
      double precision :: cc(ccl1:cch1,ccl2:cch2,ccl3:cch3)
      double precision :: ecx(ecxl1:ecxh1,ecxl2:ecxh2, ecxl3:ecxh3)
      double precision :: ecy(ecyl1:ecyh1,ecyl2:ecyh2, ecyl3:ecyh3)
      double precision :: ecz(eczl1:eczh1,eczl2:eczh2, eczl3:eczh3)
      double precision :: dx(3)

      ! Local variables
      integer          :: i,j,k

      !$OMP PARALLEL PRIVATE(i,j,k)
      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               ecx(i,j,k) = 0.5d0 * ( cc(i,j,k) - cc(i-1,j,k) ) / dx(1)
            enddo
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               ecy(i,j,k) = 0.5d0 * ( cc(i,j,k) - cc(i,j-1,k) ) / dx(2)
            enddo
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ecz(i,j,k) = 0.5d0 * ( cc(i,j,k) - cc(i,j,k-1) ) / dx(3)
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
         
      end subroutine fort_compute_ec

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
 
      subroutine ca_test_residual(lo, hi, &
           rhs, rhl1, rhl2, rhl3, rhh1, rhh2, rhh3, &
           ecx, ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3, &
           ecy, ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3, &
           ecz, eczl1, eczl2, eczl3, eczh1, eczh2, eczh3, &
           dx,problo,coord_type)

      implicit none

      integer          :: lo(3),hi(3)
      integer          :: coord_type
      integer          :: rhl1, rhl2, rhl3, rhh1, rhh2, rhh3
      integer          :: ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3
      integer          :: ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3
      integer          :: eczl1, eczl2, eczl3, eczh1, eczh2, eczh3
      double precision :: rhs(rhl1:rhh1,rhl2:rhh2,rhl3:rhh3)
      double precision :: ecx(ecxl1:ecxh1,ecxl2:ecxh2, ecxl3:ecxh3)
      double precision :: ecy(ecyl1:ecyh1,ecyl2:ecyh2, ecyl3:ecyh3)
      double precision :: ecz(eczl1:eczh1,eczl2:eczh2, eczl3:eczh3)
      double precision :: dx(3),problo(3)
     
      ! Local variables
      integer          :: i,j,k
      double precision :: lapphi

      !$OMP PARALLEL DO PRIVATE(i,j,k,lapphi) 
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               lapphi = (ecx(i+1,j,k)-ecx(i,j,k)) / dx(1) + &
                        (ecy(i,j+1,k)-ecy(i,j,k)) / dx(2) + &
                        (ecz(i,j,k+1)-ecz(i,j,k)) / dx(3)
               rhs(i,j,k) = rhs(i,j,k) - lapphi
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
 
      end subroutine ca_test_residual

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine fort_average_ec ( &
           fx, fxl1, fxl2, fxl3, fxh1, fxh2, fxh3, &
           fy, fyl1, fyl2, fyl3, fyh1, fyh2, fyh3, &
           fz, fzl1, fzl2, fzl3, fzh1, fzh2, fzh3, &
           cx, cxl1, cxl2, cxl3, cxh1, cxh2, cxh3, &
           cy, cyl1, cyl2, cyl3, cyh1, cyh2, cyh3, &
           cz, czl1, czl2, czl3, czh1, czh2, czh3, &
           lo, hi, rr)

      integer lo(3),hi(3)
      integer fxl1, fxl2, fxl3, fxh1, fxh2, fxh3
      integer fyl1, fyl2, fyl3, fyh1, fyh2, fyh3
      integer fzl1, fzl2, fzl3, fzh1, fzh2, fzh3
      integer cxl1, cxl2, cxl3, cxh1, cxh2, cxh3
      integer cyl1, cyl2, cyl3, cyh1, cyh2, cyh3
      integer czl1, czl2, czl3, czh1, czh2, czh3
      double precision fx(fxl1:fxh1,fxl2:fxh2,fxl3:fxh3)
      double precision fy(fyl1:fyh1,fyl2:fyh2,fyl3:fyh3)
      double precision fz(fzl1:fzh1,fzl2:fzh2,fzl3:fzh3)
      double precision cx(cxl1:cxh1,cxl2:cxh2,cxl3:cxh3)
      double precision cy(cyl1:cyh1,cyl2:cyh2,cyl3:cyh3)
      double precision cz(czl1:czh1,czl2:czh2,czl3:czh3)
      integer rr(3)

      integer i,j,k,n,m,facx,facy,facz

      double precision invfacyz,invfacxz,invfacxy

      facx = rr(1)
      facy = rr(2)
      facz = rr(3)

      invfacyz = 1.0d0/(facy*facz)
      invfacxz = 1.0d0/(facx*facz)
      invfacxy = 1.0d0/(facx*facy)

      !$OMP PARALLEL PRIVATE(i,j,k,m,n)
      !$OMP DO
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               cx(i,j,k) = 0.d0
               do n = 0,facy-1
                  do m = 0,facz-1
                     cx(i,j,k) = cx(i,j,k) + fx(facx*i,facy*j+n,facz*k+m)
                  end do
               end do
               cx(i,j,k) = cx(i,j,k) * invfacyz
            end do
         end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               cy(i,j,k) = 0.d0
               do n = 0,facx-1
                  do m = 0,facz-1
                     cy(i,j,k) = cy(i,j,k) + fy(facx*i+n,facy*j,facz*k+m)
                  end do
               end do
               cy(i,j,k) = cy(i,j,k) * invfacxz
            end do
         end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               cz(i,j,k) = 0.d0
               do n = 0,facx-1
                  do m = 0,facy-1
                     cz(i,j,k) = cz(i,j,k) + fz(facx*i+n,facy*j+m,facz*k)
                  end do
               end do
               cz(i,j,k) = cz(i,j,k) * invfacxy
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      end subroutine fort_average_ec

