      real function frcorr( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c The value returned by this function is the supplementary central
c   attraction required to correct the grid force for softening and/or
c   finite thickness and/or truncations and tapers.  It does not include
c   the rigid halo component which is returned separately by HALFRC.
c The value is obtained from a table which is assumed to have been
c   pre-calculated by a call to SUPSET in the start up phase.  This table
c   will not normally be changed during the simulation.
c As an abrupt truncation could result in a larger-than-required central
c   force from the active mass, the sign of this corrective force could be
c   positive, ie outwards, though it should normally be negative.
c
c The values of the input radius and the returned acceleration are in internal
c   grid units.
c This attraction declines as r^-2 outside the last tabulated radius
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
c external
      real uofgr
c
c local variables
      integer iu
      real du, u
c
      frcorr = 0
c linear interpolation in u
      if( suppl )then
        if( lsupst )then
c table independent of grid type
          u = log( 1. + r ) / alp2
          iu = u
          if( iu .lt. nrsup )then
            du = u - real( iu )
            iu = iu + 1
            frcorr = ( 1. - du ) * htab( 2, iu ) +
     +                      du   * htab( 2, iu + 1 )
          else
            frcorr = htab( 2, nrsup ) * ( htab( 1, nrsup ) / r )**2
          end if
        else
          call crash( 'FRCORR', 'Table not available' )
        end if
      end if
      return
      end
