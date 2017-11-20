      subroutine jscale( x1, x2, y1, y2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to to set the x- and y-ranges of values that correspond to
c   the selected active area of the plotting surface
c
c calling arguments
      real x1, x2, y1, y2
c
c common block
c
      include 'inc/jscmmn.f'
c
c revert (if necessary) to full plotting area
      rx1 = tx1
      rx2 = tx2
      ry1 = ty1
      ry2 = ty2
      call jsrscl( x1, x2, y1, y2 )
      return
      end
