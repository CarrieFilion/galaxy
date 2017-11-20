      subroutine orbsto( t, y )
c  Copyright (C) 2017, Jerry Sellwood
      use aarrays
      implicit none
c external subroutine used by routine ode_tab to integrate an orbit in
c   a fixed, planar potential.  It stores the values of the dependent
c   variables for use later and updates the time at which they are
c   next needed.
c
c The equations of motion are:
c         udot = frtot( r ) + Lz^2 / r^3
c         rdot = u
c         phidot = Lz / r^2
c   with y(1) = u, y(2) = r, y(3) = phi
c
c calling arguments
      real*8 t, y( 3 )
c
c common block
c
      include 'inc/orbval.f'
c
c store values
      is = is + 1
      orbtab( idim + is ) = y( 2 )
      orbtab( 2 * idim + is ) = y( 3 )
      orbtab( 3 * idim + is ) = y( 1 )
c time for next value
      t = orbtab( is + 1 )
      return
      end
