      real*8 function ajfun( m, n, r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the value of the (k,m,n) Abel-Jacobi function at the given radius
c The routine uses Kalnajs's trick (1976, ApJ v205, pp 745-750) of factoring
c   out the roots for improved precision.  The required roots and
c   coefficients are read from the appropriate input file on the first call.
c
c calling arguments
      integer m, n
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/bdpps.f'
c
c local variables
      integer i
      logical first
      real*8 r2, x
      include 'inc/pi.f'
      save first
c
      data first / .true. /
c
      if( first )then
        call ajsetup
        first = .false.
      end if
c
      r2 = r * r
c polynomial in (1 - r^2)
      x = 1.d0 - r2
      ajfun = xcoeff( m, n, 0 )
      do i = 1, basis
        ajfun = ajfun + xcoeff( m, n, i ) * x**i
      end do
c factored terms (excluding r^m)
      do i = 1, n
        ajfun = ajfun * ( r2 - root( m, n, i )**2 )
      end do
c
      if( m .gt. 0 )ajfun = r**m * ajfun
c normalise to sensible units (factors omitted from David's files)
      ajfun = pi * ajfun / sqrt( 2. )
      return
      end
