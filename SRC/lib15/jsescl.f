      subroutine jsescl( x1, x2, y1, y2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to scale the current plotting (sub-)window in such a way
c   that the number of inches or cm per unit range of x is equal to
c   the same dimension per unit range of y
c
c calling arguments
      real x1, x2, y1, y2
c
c common block
c
      include 'inc/jscmmn.f'
c
c local variable
      real s
c
      xfac = ( tx2 - tx1 ) / ( x2 - x1 )
      yfac = ( ty2 - ty1 ) / ( y2 - y1 )
c adjust frame size to make these equal
      s = min( xfac, yfac )
      if( s .eq. yfac )then
        xfac = s
        s = s * ( x2 - x1 )
        rx1 = tx1 + .5 * ( tx2 - tx1 - s )
        rx2 = rx1 + s
        ry1 = ty1
        ry2 = ty2
      else
        yfac = s
        rx1 = tx1
        rx2 = tx2
        s = s * ( y2 - y1 )
        ry1 = ty1 + .5 * ( ty2 - ty1 - s )
        ry2 = ry1 + s
      end if
      call jsrscl( x1, x2, y1, y2 )
      return
      end
