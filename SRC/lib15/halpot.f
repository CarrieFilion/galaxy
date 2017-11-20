      real function halpot( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the potential of the rigid halo component.  This is not the
c   analytic function, but a pre-tabulated function which is properly
c   defined to go to zero at infinity - possible since the halo is
c   truncated at the grid edge.  The table is created by a call to halset.
c Values outside the grid decline as 1/r
c
c calling argument
      real r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/supp.f'
c
c local variables
      integer iu
      real du, u
c
c return zero if halo is not rigid
      halpot = 0
      if( rigidh )then
c table independent of grid type
        u = log( 1. + r ) / alp2
        iu = u
        if( iu .lt. nrsup - 1 )then
c linear interpolation in u
          du = u - real( iu )
          iu = iu + 1
          halpot = ( 1. - du ) * hpot( iu ) + du * hpot( iu + 1 )
        else
c vacuum solution beyond the last measured point
          halpot = hpot( nrsup ) * htab( 1, nrsup ) / r
        end if
      end if
      return
      end
