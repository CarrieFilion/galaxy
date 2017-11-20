      subroutine tmcgrp( jst, start )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to time center coordinates (or to undo this process) for one group
c   of particles.
c There are two time centering options for leapfrog set by a logical variable:
c   if( .not. lfrv ) the positions are stored half a step ahead (APLF)
c   if( lfrv ) the velocities are stored half a step behind (RVLF)
c
c RVLF time stepping has tthe disadvantage that it requires the acceleration
c    components to be available before time centered values can be estimated
c    during analysis or created or undone (by this routine)
c
c if( start )then
c   the coordinates are set out of synch ready to integrate the model forwards
c else
c   the coordinates are restored to a synchronous state
c end if
c
c calling arguments
      integer jst
      logical start
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
c external
      logical offgrd
c
c local variables
      integer i, is
      logical off
      real fac
c
c determine sense of adjustment
      if( start )then
        fac = .5
      else
        fac = -.5
      end if
c
      if( lfrv )then
c RVLF time stepping
        do is = 1, jst
c copy positions and adjust velocities
          do i = 1, ndimen
            newc( i, is ) = oldc( i, is )
            newc( i + ndimen, is ) =
     +                       oldc( i + ndimen, is ) - fac * acc( i, is )
c initial momentum of model compensates for momentum of perturber
            if( pertbn )newc( i + ndimen, is ) =
     +              newc( i + ndimen, is ) - 2 * fac * vgal( i ) / gvfac
          end do
        end do
      else
c APLF time steping
        do is = 1, jst
c copy velocities and adjust positions
          do i = 1, ndimen
            newc( i + ndimen, is ) = oldc( i + ndimen, is )
c initial momentum of model compensates for momentum of perturber
            if( pertbn )newc( i + ndimen, is ) =
     +              newc( i + ndimen, is ) - 2 * fac * vgal( i ) / gvfac
            newc( i, is ) = oldc( i, is ) + fac * newc( i + ndimen, is )
          end do
        end do
c check for particles leaving or returning to the grid
        do is = 1, jst
          off = offgrd( is )
          if( ilist .lt. nlists )then
            if( off )call addoff( is )
          else
            label( is ) = igrd( iflag( is ) )
            if( .not. off )then
              noffm = noffm - 1
              iz( is ) = nzones
            end if
          end if
        end do
      end if
      return
      end
