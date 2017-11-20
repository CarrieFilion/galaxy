      subroutine jscoor( x, y, x1, y1 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine, redundant since v14.5, to convert input values (x,y) to
c   output values (x1,y1) that could be shifted scaled 
c   this work is now left to the PGPLOT software
c
c calling arguments
      real x, y, x1, y1
c
c$$$c common blocks
c$$$c
c$$$      include 'inc/jscmmn.f'
c$$$c
c$$$c convert to cm
c$$$      x1 = ( x - sx1 ) * xfac + rx1
c$$$      y1 = ( y - sy1 ) * yfac + ry1
c$$$c limit coordinates to current window
c$$$      x1 = max( x1, fx1 )
c$$$      x1 = min( x1, fx2 )
c$$$      y1 = max( y1, fy1 )
c$$$      y1 = min( y1, fy2 )
      x1 = x
      y1 = y
      return
      end
