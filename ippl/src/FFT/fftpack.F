c========================================================================
c 
c     The IPPL Framework
c  
c     This program was prepared by the Regents of the University of
c     Visit http://people.web.psi.ch/adelmann/ for more details
c 
c========================================================================

c=======================================================================
c============================fftpack.F==================================
c=============contains routines rffti, rfftf, rfftf,====================
c==========cffti, cfftf, cfftb, sinti, and sint (from Netlib)===========
c=======================================================================

c     Cray machines don't support double precision type.  Use real.

#if defined(IPPL_T3E)
#define FLOAT REAL
#else
#define FLOAT DOUBLE PRECISION
#endif

c---------------------------DOUBLE PRECISION----------------------------
      subroutine rffti (n,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wsave(*)
      if (n .eq. 1) return
      call rffti1 (n,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rffti1 (n,wa,xxifac)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wa(*)      ,xxifac(*)    ,ntryh(4)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      jfac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         jfac(ib+2) = jfac(ib+1)
  106 continue
      jfac(3) = 2
  107 if (nl .ne. 1) go to 104
      jfac(1) = n
      jfac(2) = nf
      tpi = 6.28318530717959
      argh = tpi/float(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = jfac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = float(ld)*argh
            fi = 0.
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue
      do 111 i = 1, 15
         xxifac(i) = rfac(i)
  111 continue
      return
      end



      subroutine rfftf (n,r,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       r(*)       ,wsave(*)
      if (n .eq. 1) return
      call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rfftf1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = jfac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call radf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call radf2 (ido,l1,ch,c,wa(iw))
         go to 110
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
  110    l2 = l1
  111 continue
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
      return
      end
      subroutine radfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
      data tpi/6.28318530717959/
      arg = tpi/float(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do 101 ik=1,idl1
         ch2(ik,1) = c2(ik,1)
  101 continue
      do 103 j=2,ip
         do 102 k=1,l1
            ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
      if (nbd .gt. l1) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108       continue
  109    continue
  110 continue
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
c
      ar1 = 1.
      ai1 = 0.
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
c
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
      go to 135
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142       continue
  143    continue
  144 continue
      return
      end
      subroutine radf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,l1,5)           ,ch(ido,5,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
      return
      end
      subroutine radf3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,3,l1)           ,cc(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
      return
      end
      subroutine radf2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,2,l1)           ,cc(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
  107 return
      end
      subroutine radf4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,l1,4)           ,ch(ido,4,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data hsqt2 /.7071067811865475/
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
  107 return
      end



      subroutine rfftb (n,r,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       r(*)       ,wsave(*)
      if (n .eq. 1) return
      call rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rfftb1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = jfac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call radb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call radb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call radb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call radb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call radbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call radbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine radbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
      data tpi/6.28318530717959/
      arg = tpi/float(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
      go to 106
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
  116 ar1 = 1.
      ai1 = 0.
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido .eq. 1) return
      do 133 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  133 continue
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
      if (nbd .gt. l1) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140       continue
  141    continue
  142 continue
  143 return
      end
      subroutine radb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
      return
      end
      subroutine radb3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
      return
      end
      subroutine radb2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
  107 return
      end
      subroutine radb4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data sqrt2 /1.414213562373095/
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
  107 return
      end



      subroutine cffti (n,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cffti1 (n,wsave(iw1),wsave(iw2))
      return
      end
      subroutine cffti1 (n,wa,xxifac)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wa(*)      ,xxifac(*)    ,ntryh(4)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      jfac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         jfac(ib+2) = jfac(ib+1)
  106 continue
      jfac(3) = 2
  107 if (nl .ne. 1) go to 104
      jfac(1) = n
      jfac(2) = nf
      tpi = 6.28318530717959
      argh = tpi/float(n)
      i = 2
      l1 = 1
      do 110 k1=1,nf
         ip = jfac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido+ido+2
         ipm = ip-1
         do 109 j=1,ipm
            i1 = i
            wa(i-1) = 1.
            wa(i) = 0.
            ld = ld+l1
            fi = 0.
            argld = float(ld)*argh
            do 108 ii=4,idot,2
               i = i+2
               fi = fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            if (ip .le. 5) go to 109
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
  109    continue
         l1 = l2
  110 continue
      do 111 i = 1, 15
         xxifac(i) = rfac(i)
  111 continue
      return
      end



      subroutine cfftf (n,c,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       c(*)       ,wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cfftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end
      subroutine cfftf1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = jfac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot
         ix3 = ix2+idot
         if (na .ne. 0) go to 101
         call passf4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call passf4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call passf2 (idot,l1,c,ch,wa(iw))
         go to 105
  104    call passf2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot
         if (na .ne. 0) go to 107
         call passf3 (idot,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call passf3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         if (na .ne. 0) go to 110
         call passf5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call passf5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call passf (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call passf (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
      n2 = n+n
      do 117 i=1,n2
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine passf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip),
     2                ch2(idl1,ip)
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 idl = 2-ido
      inc = 0
      do 116 l=2,ipph
         lc = ipp2-l
         idl = idl+ido
         do 113 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = -wa(idl)*ch2(ik,ip)
  113    continue
         idlj = idl
         inc = inc+ido
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 114 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)-wai*ch2(ik,jc)
  114       continue
  115    continue
  116 continue
      do 118 j=2,ipph
         do 117 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
      nac = 1
      if (ido .eq. 2) return
      nac = 0
      do 121 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  121 continue
      do 123 j=2,ip
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
      if (idot .gt. l1) go to 127
      idij = 0
      do 126 j=2,ip
         idij = idij+2
         do 125 i=4,ido,2
            idij = idij+2
            do 124 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  124       continue
  125    continue
  126 continue
      return
  127 idj = 2-ido
      do 130 j=2,ip
         idj = idj+ido
         do 129 k=1,l1
            idij = idj
            do 128 i=4,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  128       continue
  129    continue
  130 continue
      return
      end
      subroutine passf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,-.951056516295154,
     1-.809016994374947,-.587785252292473/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti5 = cc(2,2,k)-cc(2,5,k)
         ti2 = cc(2,2,k)+cc(2,5,k)
         ti4 = cc(2,3,k)-cc(2,4,k)
         ti3 = cc(2,3,k)+cc(2,4,k)
         tr5 = cc(1,2,k)-cc(1,5,k)
         tr2 = cc(1,2,k)+cc(1,5,k)
         tr4 = cc(1,3,k)-cc(1,4,k)
         tr3 = cc(1,3,k)+cc(1,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         ch(2,k,1) = cc(2,1,k)+ti2+ti3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,5) = cr2+ci5
         ch(2,k,2) = ci2+cr5
         ch(2,k,3) = ci3+cr4
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(2,k,4) = ci3-cr4
         ch(2,k,5) = ci2-cr5
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti5 = cc(i,2,k)-cc(i,5,k)
            ti2 = cc(i,2,k)+cc(i,5,k)
            ti4 = cc(i,3,k)-cc(i,4,k)
            ti3 = cc(i,3,k)+cc(i,4,k)
            tr5 = cc(i-1,2,k)-cc(i-1,5,k)
            tr2 = cc(i-1,2,k)+cc(i-1,5,k)
            tr4 = cc(i-1,3,k)-cc(i-1,4,k)
            tr3 = cc(i-1,3,k)+cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4+wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4-wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5+wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5-wa4(i)*dr5
  103    continue
  104 continue
      return
      end
      subroutine passf3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,-.866025403784439/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         tr2 = cc(1,2,k)+cc(1,3,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ti2 = cc(2,2,k)+cc(2,3,k)
         ci2 = cc(2,1,k)+taur*ti2
         ch(2,k,1) = cc(2,1,k)+ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
         ch(2,k,2) = ci2+cr3
         ch(2,k,3) = ci2-cr3
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            tr2 = cc(i-1,2,k)+cc(i-1,3,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,2,k)+cc(i,3,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
  103    continue
  104 continue
      return
      end
      subroutine passf2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      if (ido .gt. 2) go to 102
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
         ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
         ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
            tr2 = cc(i-1,1,k)-cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ti2 = cc(i,1,k)-cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2-wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2+wa1(i)*ti2
  103    continue
  104 continue
      return
      end
      subroutine passf4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti1 = cc(2,1,k)-cc(2,3,k)
         ti2 = cc(2,1,k)+cc(2,3,k)
         tr4 = cc(2,2,k)-cc(2,4,k)
         ti3 = cc(2,2,k)+cc(2,4,k)
         tr1 = cc(1,1,k)-cc(1,3,k)
         tr2 = cc(1,1,k)+cc(1,3,k)
         ti4 = cc(1,4,k)-cc(1,2,k)
         tr3 = cc(1,2,k)+cc(1,4,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,3) = tr2-tr3
         ch(2,k,1) = ti2+ti3
         ch(2,k,3) = ti2-ti3
         ch(1,k,2) = tr1+tr4
         ch(1,k,4) = tr1-tr4
         ch(2,k,2) = ti1+ti4
         ch(2,k,4) = ti1-ti4
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,2,k)-cc(i,4,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,4,k)-cc(i-1,2,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2+wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2-wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3+wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3-wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4+wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4-wa3(i)*cr4
  103    continue
  104 continue
      return
      end



      subroutine cfftb (n,c,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       c(*)       ,wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call cfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end
      subroutine cfftb1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = jfac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot
         ix3 = ix2+idot
         if (na .ne. 0) go to 101
         call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call passb2 (idot,l1,c,ch,wa(iw))
         go to 105
  104    call passb2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot
         if (na .ne. 0) go to 107
         call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         if (na .ne. 0) go to 110
         call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
      n2 = n+n
      do 117 i=1,n2
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine passb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip),
     2                ch2(idl1,ip)
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 idl = 2-ido
      inc = 0
      do 116 l=2,ipph
         lc = ipp2-l
         idl = idl+ido
         do 113 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = wa(idl)*ch2(ik,ip)
  113    continue
         idlj = idl
         inc = inc+ido
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 114 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+wai*ch2(ik,jc)
  114       continue
  115    continue
  116 continue
      do 118 j=2,ipph
         do 117 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
      nac = 1
      if (ido .eq. 2) return
      nac = 0
      do 121 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  121 continue
      do 123 j=2,ip
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
      if (idot .gt. l1) go to 127
      idij = 0
      do 126 j=2,ip
         idij = idij+2
         do 125 i=4,ido,2
            idij = idij+2
            do 124 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  124       continue
  125    continue
  126 continue
      return
  127 idj = 2-ido
      do 130 j=2,ip
         idj = idj+ido
         do 129 k=1,l1
            idij = idj
            do 128 i=4,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  128       continue
  129    continue
  130 continue
      return
      end
      subroutine passb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti5 = cc(2,2,k)-cc(2,5,k)
         ti2 = cc(2,2,k)+cc(2,5,k)
         ti4 = cc(2,3,k)-cc(2,4,k)
         ti3 = cc(2,3,k)+cc(2,4,k)
         tr5 = cc(1,2,k)-cc(1,5,k)
         tr2 = cc(1,2,k)+cc(1,5,k)
         tr4 = cc(1,3,k)-cc(1,4,k)
         tr3 = cc(1,3,k)+cc(1,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         ch(2,k,1) = cc(2,1,k)+ti2+ti3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,5) = cr2+ci5
         ch(2,k,2) = ci2+cr5
         ch(2,k,3) = ci3+cr4
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(2,k,4) = ci3-cr4
         ch(2,k,5) = ci2-cr5
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti5 = cc(i,2,k)-cc(i,5,k)
            ti2 = cc(i,2,k)+cc(i,5,k)
            ti4 = cc(i,3,k)-cc(i,4,k)
            ti3 = cc(i,3,k)+cc(i,4,k)
            tr5 = cc(i-1,2,k)-cc(i-1,5,k)
            tr2 = cc(i-1,2,k)+cc(i-1,5,k)
            tr4 = cc(i-1,3,k)-cc(i-1,4,k)
            tr3 = cc(i-1,3,k)+cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4-wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4+wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5-wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5+wa4(i)*dr5
  103    continue
  104 continue
      return
      end
      subroutine passb3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         tr2 = cc(1,2,k)+cc(1,3,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ti2 = cc(2,2,k)+cc(2,3,k)
         ci2 = cc(2,1,k)+taur*ti2
         ch(2,k,1) = cc(2,1,k)+ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
         ch(2,k,2) = ci2+cr3
         ch(2,k,3) = ci2-cr3
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            tr2 = cc(i-1,2,k)+cc(i-1,3,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,2,k)+cc(i,3,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
  103    continue
  104 continue
      return
      end
      subroutine passb2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      if (ido .gt. 2) go to 102
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
         ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
         ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
            tr2 = cc(i-1,1,k)-cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ti2 = cc(i,1,k)-cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2+wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2-wa1(i)*ti2
  103    continue
  104 continue
      return
      end
      subroutine passb4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti1 = cc(2,1,k)-cc(2,3,k)
         ti2 = cc(2,1,k)+cc(2,3,k)
         tr4 = cc(2,4,k)-cc(2,2,k)
         ti3 = cc(2,2,k)+cc(2,4,k)
         tr1 = cc(1,1,k)-cc(1,3,k)
         tr2 = cc(1,1,k)+cc(1,3,k)
         ti4 = cc(1,2,k)-cc(1,4,k)
         tr3 = cc(1,2,k)+cc(1,4,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,3) = tr2-tr3
         ch(2,k,1) = ti2+ti3
         ch(2,k,3) = ti2-ti3
         ch(1,k,2) = tr1+tr4
         ch(1,k,4) = tr1-tr4
         ch(2,k,2) = ti1+ti4
         ch(2,k,4) = ti1-ti4
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,4,k)-cc(i,2,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,2,k)-cc(i-1,4,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
  103    continue
  104 continue
      return
      end



      subroutine sinti (n,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wsave(*)
      data pi /3.14159265358979/
      if (n .le. 1) return
      ns2 = n/2
      np1 = n+1
      dt = pi/float(np1)
      do 101 k=1,ns2
         wsave(k) = 2.*sin(k*dt)
  101 continue
      call rffti (np1,wsave(ns2+1))
      return
      end



      subroutine sint (n,x,wsave)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       x(*)       ,wsave(*)
      np1 = n+1
      iw1 = n/2+1
      iw2 = iw1+np1
      iw3 = iw2+np1
      call sint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))
      return
      end
      subroutine sint1(n,war,was,xh,x,xxifac)
c----------------------------TJW
      implicit FLOAT (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension war(*),was(*),x(*),xh(*),xxifac(*)
      data sqrt3 /1.73205080756888/
      do 100 i=1,n
      xh(i) = war(i)
      war(i) = x(i)
  100 continue
      if (n-2) 101,102,103
  101 xh(1) = xh(1)+xh(1)
      go to 106
  102 xhold = sqrt3*(xh(1)+xh(2))
      xh(2) = sqrt3*(xh(1)-xh(2))
      xh(1) = xhold
      go to 106
  103 np1 = n+1
      ns2 = n/2
      x(1) = 0.
      do 104 k=1,ns2
         kc = np1-k
         t1 = xh(k)-xh(kc)
         t2 = was(k)*(xh(k)+xh(kc))
         x(k+1) = t1+t2
         x(kc+1) = t2-t1
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+2) = 4.*xh(ns2+1)
      call rfftf1 (np1,x,xh,war,xxifac)
      xh(1) = .5*x(1)
      do 105 i=3,n,2
         xh(i-1) = -x(i)
         xh(i) = xh(i-2)+x(i-1)
  105 continue
      if (modn .ne. 0) go to 106
      xh(n) = -x(n+1)
  106 do 107 i=1,n
      x(i) = war(i)
      war(i) = xh(i)
  107 continue
      return
      end

c---------------------------SINGLE PRECISION----------------------------
      subroutine frffti (n,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wsave(*)
      if (n .eq. 1) return
      call frffti1 (n,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine frffti1 (n,wa,xxifac)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wa(*)      ,xxifac(*)    ,ntryh(4)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      jfac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         jfac(ib+2) = jfac(ib+1)
  106 continue
      jfac(3) = 2
  107 if (nl .ne. 1) go to 104
      jfac(1) = n
      jfac(2) = nf
      tpi = 6.28318530717959
      argh = tpi/float(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = jfac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = float(ld)*argh
            fi = 0.
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue
      do 111 i = 1, 15
         xxifac(i) = rfac(i)
  111 continue
      return
      end



      subroutine frfftf (n,r,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       r(*)       ,wsave(*)
      if (n .eq. 1) return
      call frfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine frfftf1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = jfac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call fradf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call fradf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call fradf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call fradf2 (ido,l1,ch,c,wa(iw))
         go to 110
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call fradf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call fradf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call fradf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call fradf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call fradfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call fradfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
  110    l2 = l1
  111 continue
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
      return
      end
      subroutine fradfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
      data tpi/6.28318530717959/
      arg = tpi/float(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do 101 ik=1,idl1
         ch2(ik,1) = c2(ik,1)
  101 continue
      do 103 j=2,ip
         do 102 k=1,l1
            ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
      if (nbd .gt. l1) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108       continue
  109    continue
  110 continue
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
c
      ar1 = 1.
      ai1 = 0.
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
c
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
      go to 135
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142       continue
  143    continue
  144 continue
      return
      end
      subroutine fradf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,l1,5)           ,ch(ido,5,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
      return
      end
      subroutine fradf3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,3,l1)           ,cc(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
      return
      end
      subroutine fradf2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,2,l1)           ,cc(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
  107 return
      end
      subroutine fradf4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,l1,4)           ,ch(ido,4,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data hsqt2 /.7071067811865475/
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
  107 return
      end



      subroutine frfftb (n,r,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       r(*)       ,wsave(*)
      if (n .eq. 1) return
      call frfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine frfftb1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = jfac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call fradb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call fradb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call fradb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call fradb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call fradb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call fradb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call fradb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call fradb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call fradbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call fradbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine fradbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
      data tpi/6.28318530717959/
      arg = tpi/float(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
      go to 106
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
  116 ar1 = 1.
      ai1 = 0.
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido .eq. 1) return
      do 133 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  133 continue
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
      if (nbd .gt. l1) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140       continue
  141    continue
  142 continue
  143 return
      end
      subroutine fradb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
      return
      end
      subroutine fradb3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
      return
      end
      subroutine fradb2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
  107 return
      end
      subroutine fradb4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data sqrt2 /1.414213562373095/
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
  107 return
      end



      subroutine fcffti (n,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call fcffti1 (n,wsave(iw1),wsave(iw2))
      return
      end
      subroutine fcffti1 (n,wa,xxifac)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wa(*)      ,xxifac(*)    ,ntryh(4)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      jfac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         jfac(ib+2) = jfac(ib+1)
  106 continue
      jfac(3) = 2
  107 if (nl .ne. 1) go to 104
      jfac(1) = n
      jfac(2) = nf
      tpi = 6.28318530717959
      argh = tpi/float(n)
      i = 2
      l1 = 1
      do 110 k1=1,nf
         ip = jfac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido+ido+2
         ipm = ip-1
         do 109 j=1,ipm
            i1 = i
            wa(i-1) = 1.
            wa(i) = 0.
            ld = ld+l1
            fi = 0.
            argld = float(ld)*argh
            do 108 ii=4,idot,2
               i = i+2
               fi = fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            if (ip .le. 5) go to 109
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
  109    continue
         l1 = l2
  110 continue
      do 111 i = 1, 15
         xxifac(i) = rfac(i)
  111 continue
      return
      end



      subroutine fcfftf (n,c,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       c(*)       ,wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call fcfftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end
      subroutine fcfftf1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = jfac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot
         ix3 = ix2+idot
         if (na .ne. 0) go to 101
         call fpassf4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call fpassf4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call fpassf2 (idot,l1,c,ch,wa(iw))
         go to 105
  104    call fpassf2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot
         if (na .ne. 0) go to 107
         call fpassf3 (idot,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call fpassf3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         if (na .ne. 0) go to 110
         call fpassf5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call fpassf5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call fpassf (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call fpassf (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
      n2 = n+n
      do 117 i=1,n2
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine fpassf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip),
     2                ch2(idl1,ip)
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 idl = 2-ido
      inc = 0
      do 116 l=2,ipph
         lc = ipp2-l
         idl = idl+ido
         do 113 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = -wa(idl)*ch2(ik,ip)
  113    continue
         idlj = idl
         inc = inc+ido
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 114 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)-wai*ch2(ik,jc)
  114       continue
  115    continue
  116 continue
      do 118 j=2,ipph
         do 117 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
      nac = 1
      if (ido .eq. 2) return
      nac = 0
      do 121 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  121 continue
      do 123 j=2,ip
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
      if (idot .gt. l1) go to 127
      idij = 0
      do 126 j=2,ip
         idij = idij+2
         do 125 i=4,ido,2
            idij = idij+2
            do 124 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  124       continue
  125    continue
  126 continue
      return
  127 idj = 2-ido
      do 130 j=2,ip
         idj = idj+ido
         do 129 k=1,l1
            idij = idj
            do 128 i=4,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  128       continue
  129    continue
  130 continue
      return
      end
      subroutine fpassf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,-.951056516295154,
     1-.809016994374947,-.587785252292473/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti5 = cc(2,2,k)-cc(2,5,k)
         ti2 = cc(2,2,k)+cc(2,5,k)
         ti4 = cc(2,3,k)-cc(2,4,k)
         ti3 = cc(2,3,k)+cc(2,4,k)
         tr5 = cc(1,2,k)-cc(1,5,k)
         tr2 = cc(1,2,k)+cc(1,5,k)
         tr4 = cc(1,3,k)-cc(1,4,k)
         tr3 = cc(1,3,k)+cc(1,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         ch(2,k,1) = cc(2,1,k)+ti2+ti3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,5) = cr2+ci5
         ch(2,k,2) = ci2+cr5
         ch(2,k,3) = ci3+cr4
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(2,k,4) = ci3-cr4
         ch(2,k,5) = ci2-cr5
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti5 = cc(i,2,k)-cc(i,5,k)
            ti2 = cc(i,2,k)+cc(i,5,k)
            ti4 = cc(i,3,k)-cc(i,4,k)
            ti3 = cc(i,3,k)+cc(i,4,k)
            tr5 = cc(i-1,2,k)-cc(i-1,5,k)
            tr2 = cc(i-1,2,k)+cc(i-1,5,k)
            tr4 = cc(i-1,3,k)-cc(i-1,4,k)
            tr3 = cc(i-1,3,k)+cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4+wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4-wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5+wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5-wa4(i)*dr5
  103    continue
  104 continue
      return
      end
      subroutine fpassf3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,-.866025403784439/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         tr2 = cc(1,2,k)+cc(1,3,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ti2 = cc(2,2,k)+cc(2,3,k)
         ci2 = cc(2,1,k)+taur*ti2
         ch(2,k,1) = cc(2,1,k)+ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
         ch(2,k,2) = ci2+cr3
         ch(2,k,3) = ci2-cr3
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            tr2 = cc(i-1,2,k)+cc(i-1,3,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,2,k)+cc(i,3,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
  103    continue
  104 continue
      return
      end
      subroutine fpassf2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      if (ido .gt. 2) go to 102
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
         ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
         ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
            tr2 = cc(i-1,1,k)-cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ti2 = cc(i,1,k)-cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2-wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2+wa1(i)*ti2
  103    continue
  104 continue
      return
      end
      subroutine fpassf4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti1 = cc(2,1,k)-cc(2,3,k)
         ti2 = cc(2,1,k)+cc(2,3,k)
         tr4 = cc(2,2,k)-cc(2,4,k)
         ti3 = cc(2,2,k)+cc(2,4,k)
         tr1 = cc(1,1,k)-cc(1,3,k)
         tr2 = cc(1,1,k)+cc(1,3,k)
         ti4 = cc(1,4,k)-cc(1,2,k)
         tr3 = cc(1,2,k)+cc(1,4,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,3) = tr2-tr3
         ch(2,k,1) = ti2+ti3
         ch(2,k,3) = ti2-ti3
         ch(1,k,2) = tr1+tr4
         ch(1,k,4) = tr1-tr4
         ch(2,k,2) = ti1+ti4
         ch(2,k,4) = ti1-ti4
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,2,k)-cc(i,4,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,4,k)-cc(i-1,2,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2+wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2-wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3+wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3-wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4+wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4-wa3(i)*cr4
  103    continue
  104 continue
      return
      end



      subroutine fcfftb (n,c,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       c(*)       ,wsave(*)
      if (n .eq. 1) return
      iw1 = n+n+1
      iw2 = iw1+n+n
      call fcfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
      return
      end
      subroutine fcfftb1 (n,c,ch,wa,xxifac)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(*)      ,c(*)       ,wa(*)      ,xxifac(*)
      dimension rfac(15), jfac(15)
      equivalence (rfac(1), jfac(1))
      do 100 i = 1, 15
         rfac(i) = xxifac(i)
  100 continue
      nf = jfac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = jfac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot
         ix3 = ix2+idot
         if (na .ne. 0) go to 101
         call fpassb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call fpassb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call fpassb2 (idot,l1,c,ch,wa(iw))
         go to 105
  104    call fpassb2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot
         if (na .ne. 0) go to 107
         call fpassb3 (idot,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call fpassb3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         if (na .ne. 0) go to 110
         call fpassb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call fpassb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call fpassb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call fpassb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
      n2 = n+n
      do 117 i=1,n2
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine fpassb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip),
     2                ch2(idl1,ip)
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 idl = 2-ido
      inc = 0
      do 116 l=2,ipph
         lc = ipp2-l
         idl = idl+ido
         do 113 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = wa(idl)*ch2(ik,ip)
  113    continue
         idlj = idl
         inc = inc+ido
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 114 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+wai*ch2(ik,jc)
  114       continue
  115    continue
  116 continue
      do 118 j=2,ipph
         do 117 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
      nac = 1
      if (ido .eq. 2) return
      nac = 0
      do 121 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  121 continue
      do 123 j=2,ip
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
      if (idot .gt. l1) go to 127
      idij = 0
      do 126 j=2,ip
         idij = idij+2
         do 125 i=4,ido,2
            idij = idij+2
            do 124 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  124       continue
  125    continue
  126 continue
      return
  127 idj = 2-ido
      do 130 j=2,ip
         idj = idj+ido
         do 129 k=1,l1
            idij = idj
            do 128 i=4,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  128       continue
  129    continue
  130 continue
      return
      end
      subroutine fpassb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti5 = cc(2,2,k)-cc(2,5,k)
         ti2 = cc(2,2,k)+cc(2,5,k)
         ti4 = cc(2,3,k)-cc(2,4,k)
         ti3 = cc(2,3,k)+cc(2,4,k)
         tr5 = cc(1,2,k)-cc(1,5,k)
         tr2 = cc(1,2,k)+cc(1,5,k)
         tr4 = cc(1,3,k)-cc(1,4,k)
         tr3 = cc(1,3,k)+cc(1,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         ch(2,k,1) = cc(2,1,k)+ti2+ti3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,5) = cr2+ci5
         ch(2,k,2) = ci2+cr5
         ch(2,k,3) = ci3+cr4
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(2,k,4) = ci3-cr4
         ch(2,k,5) = ci2-cr5
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti5 = cc(i,2,k)-cc(i,5,k)
            ti2 = cc(i,2,k)+cc(i,5,k)
            ti4 = cc(i,3,k)-cc(i,4,k)
            ti3 = cc(i,3,k)+cc(i,4,k)
            tr5 = cc(i-1,2,k)-cc(i-1,5,k)
            tr2 = cc(i-1,2,k)+cc(i-1,5,k)
            tr4 = cc(i-1,3,k)-cc(i-1,4,k)
            tr3 = cc(i-1,3,k)+cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4-wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4+wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5-wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5+wa4(i)*dr5
  103    continue
  104 continue
      return
      end
      subroutine fpassb3 (ido,l1,cc,ch,wa1,wa2)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         tr2 = cc(1,2,k)+cc(1,3,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ti2 = cc(2,2,k)+cc(2,3,k)
         ci2 = cc(2,1,k)+taur*ti2
         ch(2,k,1) = cc(2,1,k)+ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
         ch(2,k,2) = ci2+cr3
         ch(2,k,3) = ci2-cr3
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            tr2 = cc(i-1,2,k)+cc(i-1,3,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,2,k)+cc(i,3,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
  103    continue
  104 continue
      return
      end
      subroutine fpassb2 (ido,l1,cc,ch,wa1)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      if (ido .gt. 2) go to 102
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
         ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
         ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
            tr2 = cc(i-1,1,k)-cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ti2 = cc(i,1,k)-cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2+wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2-wa1(i)*ti2
  103    continue
  104 continue
      return
      end
      subroutine fpassb4 (ido,l1,cc,ch,wa1,wa2,wa3)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti1 = cc(2,1,k)-cc(2,3,k)
         ti2 = cc(2,1,k)+cc(2,3,k)
         tr4 = cc(2,4,k)-cc(2,2,k)
         ti3 = cc(2,2,k)+cc(2,4,k)
         tr1 = cc(1,1,k)-cc(1,3,k)
         tr2 = cc(1,1,k)+cc(1,3,k)
         ti4 = cc(1,2,k)-cc(1,4,k)
         tr3 = cc(1,2,k)+cc(1,4,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,3) = tr2-tr3
         ch(2,k,1) = ti2+ti3
         ch(2,k,3) = ti2-ti3
         ch(1,k,2) = tr1+tr4
         ch(1,k,4) = tr1-tr4
         ch(2,k,2) = ti1+ti4
         ch(2,k,4) = ti1-ti4
  101 continue
      return
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,4,k)-cc(i,2,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,2,k)-cc(i-1,4,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
  103    continue
  104 continue
      return
      end



      subroutine fsinti (n,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       wsave(*)
      data pi /3.14159265358979/
      if (n .le. 1) return
      ns2 = n/2
      np1 = n+1
      dt = pi/float(np1)
      do 101 k=1,ns2
         wsave(k) = 2.*sin(k*dt)
  101 continue
      call frffti (np1,wsave(ns2+1))
      return
      end



      subroutine fsint (n,x,wsave)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension       x(*)       ,wsave(*)
      np1 = n+1
      iw1 = n/2+1
      iw2 = iw1+np1
      iw3 = iw2+np1
      call fsint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))
      return
      end
      subroutine fsint1(n,war,was,xh,x,xxifac)
c----------------------------TJW
      implicit real (a-h,o-z)
      implicit integer (i-n)
c----------------------------TJW
      dimension war(*),was(*),x(*),xh(*),xxifac(*)
      data sqrt3 /1.73205080756888/
      do 100 i=1,n
      xh(i) = war(i)
      war(i) = x(i)
  100 continue
      if (n-2) 101,102,103
  101 xh(1) = xh(1)+xh(1)
      go to 106
  102 xhold = sqrt3*(xh(1)+xh(2))
      xh(2) = sqrt3*(xh(1)-xh(2))
      xh(1) = xhold
      go to 106
  103 np1 = n+1
      ns2 = n/2
      x(1) = 0.
      do 104 k=1,ns2
         kc = np1-k
         t1 = xh(k)-xh(kc)
         t2 = was(k)*(xh(k)+xh(kc))
         x(k+1) = t1+t2
         x(kc+1) = t2-t1
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+2) = 4.*xh(ns2+1)
      call frfftf1 (np1,x,xh,war,xxifac)
      xh(1) = .5*x(1)
      do 105 i=3,n,2
         xh(i-1) = -x(i)
         xh(i) = xh(i-2)+x(i-1)
  105 continue
      if (modn .ne. 0) go to 106
      xh(n) = -x(n+1)
  106 do 107 i=1,n
      x(i) = war(i)
      war(i) = xh(i)
  107 continue
      return
      end



c========================================================================
c $RCSfile: fftpack.F,v $   $Author: adelmann $
c $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
c IPPL_VERSION_ID: $Id: fftpack.F,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
c========================================================================

