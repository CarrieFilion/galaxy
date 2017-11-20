      real*8 function sphrho( x, y, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Calculates the volume density at an arbitrary point by integrating the
c   DF over all allowed velocities.  The position coordinates and
c   returned value are in natural units.
c
c The routine branches depending upon whether the DF depends only on E
c   or both E and h
c
c calling arguments in natural units
      real*8 x, y, z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
c externals
      real*8 sphrfn, phitot, phisph, velif1
c
c local arrays
      real*8, allocatable :: absc( : ), wt( : )
c
c local variables
      integer i, n
      real*8 r, umax, zero
      parameter ( n = 16, zero = 0 )
      include 'inc/pi.f'
c
c get potential
      if( ( cdft( icmp ) .eq. 'DJDZ' ) .or.
     +    ( cdft( icmp ) .eq. 'DFIT' ) .or.
     +    ( cdft( icmp ) .eq. 'ADIA' ) )then
        if( sphrod( icmp ) )then
          r = sqrt( x * x + y * y )
          pot = phisph( r, z )
        else
          r = sqrt( x * x + y * y + z * z )
          pot = phisph( r, zero )
        end if
      else
        r = sqrt( x * x + y * y + z * z )
        pot = phitot( r )
      end if
c check that some bound orbits exist
      sphrho = 0.
      if( pot .lt. phimax )then
c integrate over velocities
        umax = sqrt( 2. * ( phimax - pot ) )
        if( cdft( icmp ) .eq. 'DJDZ' )then
c f(E,h) - therefore 2-D quadrature
          iu = 0
          iv = 0
          rad = r
          allocate ( absc( n ) )
          allocate ( wt( n ) )
c non-adaptive Gauss-Legendre quadrature for outer integral
          call GLqtab( zero, umax, n, wt, absc )
          do i = 1, n
            r = absc( i )
            sphrho = sphrho + wt( i ) * velif1( r )
          end do
        else
c f(E) only - therefore only 1-D quadrature
          allocate ( absc( n ) )
          allocate ( wt( n ) )
          call GLqtab( zero, umax, n, wt, absc )
          do i = 1, n
            r = absc( i )
            sphrho = sphrho + wt( i ) * sphrfn( r )
          end do
          sphrho = 4 * pi * sphrho
        end if
      end if
      return
      end
