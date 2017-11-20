      subroutine polcat
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to convert acceleration components in polar coordinates
c   to Cartesian components at every grid point
      use aarrays
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
      integer i, j, l, n
      real cp, fx, fy, sp, t
      include 'inc/pi.f'
c
      if( p2d .or. p3d )then
        n = nr( jgrid )
        if( p3d )n = nr( jgrid ) * ngz
c outer loop over spokes because angles are constant
        do j = 1, na
          t = 2. * pi * real( j - 1 ) / real( na )
          cp = cos( t )
          sp = sin( t )
          l = j
c work over rings (and planes)
          do i = 1, n
            fx = grdfld( l, 1 ) * cp - grdfld( l, 2 ) * sp
            fy = grdfld( l, 1 ) * sp + grdfld( l, 2 ) * cp
            grdfld( l, 1 ) = fx
            grdfld( l, 2 ) = fy
c update pointers
            l = l + na
          end do
c end loop over spokes
        end do
      else
        call crash( 'POLCAT', 'Unrecognized grid' )
      end if
      return
      end
