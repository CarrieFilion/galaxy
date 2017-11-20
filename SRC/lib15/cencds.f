      subroutine cencds( jst, coords )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c makes a first-order estimate of the time-centered coordinates of all
c    particles in the current group
c For RVLF time centering in leap-frog, the time-centered velocity is
c    advanced half a step
c For APLF time centering, the position is backed up half a step
c
c Called from ANLGRP - which is called right after the group has been
c   accelerated.  This essential for RVLF time centering because that
c   is the only moment when the velocities at both times are available
c
c calling arguments
      integer jst
      real coords( 6, jst )
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
      real tsfac
c
c local variables
      integer i, is
      real tfac
c
      if( lfrv )then
c synchronize velocities
        do is = 1, jst
c get time-step factor
          tfac = tsfac( iz( is ) )
c compute time-centered components
          do i = 1, ndimen
            coords( i, is ) = oldc( i, is )
            coords( i + ndimen, is ) =
     +             ( oldc( i + ndimen, is ) + .5 * acc( i, is ) ) / tfac
          end do
        end do
      else
c synchronize positions
        do is = 1, jst
c get time-step factor
          tfac = tsfac( iz( is ) )
c compute time-centered components
          do i = 1, ndimen
            coords( i, is ) =
     +                       oldc( i, is ) - .5 * oldc( i + ndimen, is )
            coords( i + ndimen, is ) = oldc( i + ndimen, is ) / tfac
          end do
        end do
      end if
      return
      end
