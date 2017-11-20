      subroutine jspnta( x, y, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to output a point of the currently requred size at the
c   point (x,y)
c   the argument n is no longer used
c
c calling arguments
      integer n
      real x, y
c
c common blocks
c
      include 'inc/jscmmn.f'
c
      real xs, ys
c
c convert to device units
      xs = x
      ys = y
c pgplot
      call pgpoint( 1, xs, ys, -1 )
      return
      end
