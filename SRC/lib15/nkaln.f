      real*8 function NKaln( m, alpha )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c kernel for solution of potential by logarithmic spiral expansion.
c   It returns the reduced potential of a unit amplitude log spiral
c   source distribution evaluated in the plane z=0 as given by formula
c   (12) of Kalnajs (1971 ApJ v166, p275), but multiplied by 2pi.  The
c   corresponding formula (2-188) of Binney & Tremaine (p83) is WRONG.
c Uses formula 6.1.25 from Abramowitz & Stegun for the product of two
c   Gamma functions having complex conjugate arguments
c
c calling arguments
      integer m
      real*8 alpha
c
c external
      real*8 Gammaf
c
c local variables
      integer n
      real*8 part, ratio, term, xden, xnum, y2
      include 'inc/pi.f'
c
      if( m .lt. 0 )call crash( 'NKALN', 'Invalid calling argument' )
c real parts of Gamma function arguments
      xnum = .5 * ( dble( m ) + .5 )
      xden = .5 * ( dble( m ) + 1.5 )
      y2 = ( .5 * alpha )**2
c infinite product part of solution
      term = ( 1. + y2 / xden**2 ) / ( 1. + y2 / xnum**2 )
      part = term
      n = 0
      do while ( abs( term - 1. ) .gt. 1.d-8 )
        n = n + 1
        term = ( 1. + y2 / ( xden + dble( n ) )**2 ) /
     +         ( 1. + y2 / ( xnum + dble( n ) )**2 )
        part = part * term
      end do
c ratio of Gamma functions part: [ Gamma( xnum ) / Gamma( xden ) ]**2
      ratio = 1
      do while ( xnum .gt. 1.d0 )
        xnum = xnum - 1
        xden = xden - 1
        ratio = ratio * ( xnum / xden )
      end do
      n = 0
      ratio = ratio * Gammaf( xnum, n ) / Gammaf( xden, n )
      ratio = ratio**2
c final result
      NKaln = pi * ratio * part
      return
      end
