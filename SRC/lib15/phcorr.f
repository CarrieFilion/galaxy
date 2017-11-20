      real function phcorr( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c The value returned by this function is the potential whose derivative
c   yields FRCORR.  It does not include the rigid halo contribution which
c   is returned separately by HALPOT
c The value is obtained from a table which is assumed to have been
c   pre-calculated by a call to SUPSET in the start up phase and the table
c   will not normally change during the simulation
c
c The values of the input radius and the returned potential are in internal
c   grid units.
c The potential declines as r^-1 outside the last tabulated radius
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
c this piece of garbage code is inserted to defeat the optimizer, which
c   otherwise causes a floating point exception in this routine
      integer icall
      save icall
      data icall / 0 /
      icall = icall + 1
      icall = mod( icall, 1000 )
      if( icall .lt. 0 )print *, 'starting phcorr', r
c
c return zero if no supplementary forces
      phcorr = 0
      if( suppl )then
        if( lsupst )then
c table independent of grid type
          u = log( 1. + r ) / alp2
c linear interpolation in u
          iu = u
          if( iu .lt. nrsup )then
            du = u - real( iu )
            iu = iu + 1
            phcorr = ( 1. - du ) * htab( 3, iu ) +
     +                      du   * htab( 3, iu + 1 )
          else
            phcorr = htab( 3, nrsup ) * htab( 1, nrsup ) / r
          end if
        else
          call crash( 'PHCORR', 'Table not available' )
        end if
      end if
c      if( mod( icall, 1000 ) .eq. 0 )print *, '    done phcorr', phcorr
      return
      end
