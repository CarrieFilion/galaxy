      subroutine jssize( x1, x2, y1, y2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to specify the rectangular fractional part of the current
c   plotting device surface that will be enclosed by the
c   top and bottom x-axes and the left and right y-axes
c
c calling arguments
      real x1, x2, y1, y2
c
c common blocks
c
      include 'inc/jscmmn.f'
c
      tx1 = x1 * fx2 + ( 1. - x1 ) * fx1
      tx2 = x2 * fx2 + ( 1. - x2 ) * fx1
      ty1 = y1 * fy2 + ( 1. - y1 ) * fy1
      ty2 = y2 * fy2 + ( 1. - y2 ) * fy1
      rx1 = tx1
      rx2 = tx2
      ry1 = ty1
      ry2 = ty2
c pgplot dimensions are inches
      call pgvsiz( rx1 / 2.54, rx2 / 2.54, ry1 / 2.54, ry2 / 2.54 )
      return
      end
