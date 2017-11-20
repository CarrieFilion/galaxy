      real*8 function Gam1o2( x, y, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c evaluates the ratio Gamma(x)/[Gamma(y)Gamma(z)}
c
c calling arguments
      real*8 x, y, z
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
      integer i, ifail, l, m, n
      real*8 a, b, c, fac
c
c convert to factorials
      a = x - 1.
      b = max( y, z ) - 1.
      c = min( y, z ) - 1.
c find highest integers
      l = a
      m = b
      n = c
c adjust for numerators smaller than denominators
      if( n .gt. l )then
        n = 0
        m = l
      else if( n .gt. l - m )then
        n = l - m
      end if
      n = max( n, 0 )
      m = max( m, 0 )
c irreducible part
      fac = a - dble( m + n - 1 )
      Gam1o2 = Gammaf( fac, ifail )
      fac = b - dble( m - 1 )
      Gam1o2 = Gam1o2 / Gammaf( fac, ifail )
      fac = c - dble( n - 1 )
      Gam1o2 = Gam1o2 / Gammaf( fac, ifail )
c factorials
      if( m .gt. 0 )then
        do i = 1, m
          Gam1o2 = Gam1o2 * ( a / b )
          a = a - 1.
          b = b - 1.
        end do
      end if
      if( n .gt. 0 )then
        do i = 1, n
          Gam1o2 = Gam1o2 * ( a / c )
          a = a - 1.
          c = c - 1.
        end do
      end if
      return
      end
