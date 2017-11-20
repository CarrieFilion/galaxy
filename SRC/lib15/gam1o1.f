      real*8 function Gam1o1( x, y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c evaluates the ratio Gamma(x)/Gamma(y)
c
c calling arguments
      real*8 x, y
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      real*8 Gammaf
c
c local variables
      integer i, ifail, n
      real*8 a, b, fac
c
c convert to factorials
      a = x - 1.
      b = y - 1.
c find highest integer that can be factored out
      n = min( a, b )
      n = max( n, 0 )
c irreducible part
      fac = x - dble( n )
      Gam1o1 = Gammaf( fac, ifail )
      fac = y - dble( n )
      Gam1o1 = Gam1o1 / Gammaf( fac, ifail )
c factorials
      if( n .gt. 0 )then
        do i = 1, n
          Gam1o1 = Gam1o1 * ( a / b )
          a = a - 1.
          b = b - 1.
        end do
      end if
      return
      end
