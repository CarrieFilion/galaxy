      real*8 function ajder( m, n, r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the radial derivative of the (k,m,n) Abel-Jacobi function at the
c   given radius
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
c external
      real*8 ajfun
c
c local variables
      integer i, j, k, jj
      logical first
      real*8 r2, x, tmp, xpoly, xpolyder, facpart, facpartder,
     +       rm, rmder
      save first
      include 'inc/pi.f'
c
      data first / .true. /
c
      if( first )then
        tmp = ajfun( m, n, r )
        first = .false.
      end if

      k = basis
      r2 = r*r
      x = 1. - r2

* Compute r^m and its derivative

      if( m .eq. 0 )then
        rm = 1
        rmder = 0.
      else
        rm = r**m
        rmder = real( m ) * r**( m - 1 )
      end if

* Compute factored part of ajfun (not including r^m)

      facpart = 1
      do j = 1, n
        facpart = facpart * ( r2 - root( m, n, j )**2 )
      end do

* Compute derivative of factored part

      facpartder = 0.
      do j = 1, n
        tmp = 1
        do jj = 1, j-1
          tmp = tmp * ( r2 - root( m, n, jj )**2 )
        end do
        tmp = 2. * r * tmp
        do jj = j+1, n
          tmp = tmp * ( r2 - root( m, n, jj )**2 )
        end do
        facpartder = facpartder + tmp
      end do

* Compute polynomial in (1 - r^2)

      xpoly = xcoeff( m, n, 0 )
      do i = 1, k
        xpoly = xpoly + xcoeff( m, n, i ) * x**i
      end do

* Compute derivative of polynomial in (1 - r^2)

      if( k .eq. 0 )then
        xpolyder = 0.
      else
        xpolyder = xcoeff( m, n, 1 )
        do i = 2, k
          xpolyder = xpolyder +
     &            real( i ) * xcoeff( m, n, i ) * x**( i - 1 )
        end do
        xpolyder = -2. * r * xpolyder
      end if

* Finally calculate full required derivative

      ajder = rmder * facpart * xpoly
     &       + rm * facpartder * xpoly
     &       + rm * facpart * xpolyder

c normalise to sensible units (factors omitted from David's files)
      ajder = pi * ajder / sqrt( 2. )

      return
      end
