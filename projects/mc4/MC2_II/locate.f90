      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n,i
      REAL x,xx(n)
      INTEGER jl,jm,ju

!hpf$ distribute (block) :: xx

      jl=0
      ju=n+1

10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
