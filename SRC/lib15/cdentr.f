      subroutine cdentr( t, u, x, y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to transform coordinates from the polar array used by drmode
c   to Cartesians for plotting
c
c calling arguments
      real t, u, x, y
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      real rad( 1 )
      common / trcdent / rad
c
c local variables
      integer i
      real r, th
c
      th = alpha * ( t - 1. )
      i = u + uoffst
      x = u + uoffst - real( i )
      r = ( 1. - x ) * rad( i ) + x * rad( i + 1 )
c radius should be greater than 0
      if( ( u .gt. 1. ) .and. ( r .le. 0. ) )then
        print *, u, i, rad( i ), rad( i + 1 ), r
        call crash( 'CDENTR', 'Radius array not set?' )
      end if
      x = r * cos( th )
      y = r * sin( th )
      return
      end
