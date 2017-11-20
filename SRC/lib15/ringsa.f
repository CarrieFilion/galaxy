      subroutine ringsa( jst, coords, nwring, wring )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Stores positions from any of current group of particles that are members
c   of the rings of test particles
c Called from ANLGRP
c
c calling arguments
      integer jst, nwring
      real coords( 6, jst ), wring( nwring )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/buffer.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, is, j, k, l
c
c base address of ring particles
      j = nwpp * ( nbod - nsp( ncmp ) )
c work through group
      do is = 1, jst
        if( iflag( is ) .eq. ncmp )then
          k = ndimen * ( ( loc( is ) - j ) / nwpp )
c check space
          if( ( k .lt. 0 ) .or. ( k + l .gt. nwring ) )then
            print *, is, loc( is ), j, k, l
            call crash( 'RINGSA', 'Out of array bound' )
          end if
c
          do i = 1, l
            wring( k + i ) = coords( i, is )
          end do
        end if
      end do
      return
      end
