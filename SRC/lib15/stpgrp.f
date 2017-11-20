      subroutine stpgrp( jst )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c advances the motion of the current group of particles
c
c calling argument
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c external
      logical offgrd
c
c local variables
      integer i, is
      logical off
c
c time-centred leap-frog
      do is = 1, jst
        do i = 1, ndimen
          newc( i + ndimen, is ) = oldc( i + ndimen, is ) + acc( i, is )
          newc( i, is ) = oldc( i, is ) + newc( i + ndimen, is )
        end do
      end do
c shift grid relative to particles if required
      if( gshift )then
        if( izone .ne. 1 )call crash( 'STPGRP',
     +                 'Shifts for multiple time steps not programmed' )
        do is = 1, jst
          do i = 1, ndimen
            newc( i, is ) = newc( i, is ) + xshift( i )
          end do
        end do
      end if
c check for particles leaving or returning to the grid
      do is = 1, jst
        off = .false.
        if( .not. stphev )off = offgrd( is )
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
      return
      end
