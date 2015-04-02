c========================================================================
c 
c     The IPPL Framework
c  
c     This program was prepared by the Regents of the University of
c     Visit http://people.web.psi.ch/adelmann/ for more details
c 
c========================================================================

c========================================================================
c=========================SCSL_FFT.F=====================================
c========================================================================
c     I'm adding some Fortran wrappers around SGI/Cray Scientific Library
c     routines, so we can call them from C++ code
c========================================================================


c     The Cray T3E does not have double-precision versions of these 
c     routines, since its routines are double-precision by default.

#ifndef IPPL_T3E

      subroutine sgizzfft(isign,n,scale,x,y,table,work,isys)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c------------------------------------------------------------------------
      dimension x(*), y(*), table(*), work(*)
      call zzfft(isign,n,scale,x,y,table,work,isys)
      return
      end

      subroutine sgidzfft(isign,n,scale,x,y,table,work,isys)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c------------------------------------------------------------------------
      dimension x(*), y(*), table(*), work(*)
      call dzfft(isign,n,scale,x,y,table,work,isys)
      return
      end

      subroutine sgizdfft(isign,n,scale,x,y,table,work,isys)
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c------------------------------------------------------------------------
      dimension x(*), y(*), table(*), work(*)
      call zdfft(isign,n,scale,x,y,table,work,isys)
      return
      end

#endif


c     Single-precision routines

      subroutine sgiccfft(isign,n,scale,x,y,table,work,isys)
c------------------------------------------------------------------------
      implicit real (a-h,o-z)
      implicit integer (i-n)
c------------------------------------------------------------------------
      dimension x(*), y(*), table(*), work(*)
      call ccfft(isign,n,scale,x,y,table,work,isys)
      return
      end

      subroutine sgiscfft(isign,n,scale,x,y,table,work,isys)
c------------------------------------------------------------------------
      implicit real (a-h,o-z)
      implicit integer (i-n)
c------------------------------------------------------------------------
      dimension x(*), y(*), table(*), work(*)
      call scfft(isign,n,scale,x,y,table,work,isys)
      return
      end

      subroutine sgicsfft(isign,n,scale,x,y,table,work,isys)
c------------------------------------------------------------------------
      implicit real (a-h,o-z)
      implicit integer (i-n)
c------------------------------------------------------------------------
      dimension x(*), y(*), table(*), work(*)
      call csfft(isign,n,scale,x,y,table,work,isys)
      return
      end


c========================================================================
c $RCSfile: SCSL_FFT.F,v $   $Author: adelmann $
c $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
c IPPL_VERSION_ID: $Id: SCSL_FFT.F,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
c========================================================================
