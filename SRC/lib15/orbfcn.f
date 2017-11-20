      subroutine orbfcn( t, y, f )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrator for an orbit in a fixed, planar potential
c   external subroutine used by NAg routine D02CJF (Adams)
c   The equations of motion are:
c         udot = frtot( r ) + Lz^2 / r^3
c         rdot = u
c         phidot = Lz / r^2
c   with y(1) = u, y(2) = r, y(3) = phi
c
c calling arguments
      real*8 t, y( 3 ), f( 3 )
c
c common block
c
      include 'inc/orbval.f'
c
c external
      real*8 frtot
c
      f( 1 ) = frtot( y( 2 ) ) + crrntL**2 / y( 2 )**3
      f( 2 ) = y( 1 )
      f( 3 ) = crrntL / ( y( 2 ) * y( 2 ) )
      return
      end
