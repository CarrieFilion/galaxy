      subroutine bhboun
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
c computes the bounding values of the smallest Cartesian box that will
c   enclose all the particles.  Needed for the Barnes-Hut tree algrithm
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variables
      integer i, ip, jp
c
c initialise range bounds
      do i = 1, ndimen
        bound( 1, i ) = 1000
        bound( 2, i ) = -1000
      end do
c work through particles to find the extent of the occupied region
      do ip = 1, lpf, nwpp
        do i = 1, ndimen
          jp = ip - 1 + i
          bound( 1, i ) = min( bound( 1, i ), ptcls( jp ) )
          bound( 2, i ) = max( bound( 2, i ), ptcls( jp ) )
        end do
      end do
      return
      end
