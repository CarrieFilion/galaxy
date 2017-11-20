      subroutine jsarrw( x1, y1, x2, y2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to draw an arrow from the point (x1,y1) to (x2,y2)
c
c calling arguments
      real x1, y1, x2, y2
c
c common block
c
      include 'inc/jscmmn.f'
c
c local variables
      real dx, dy, headl, x, y
      parameter ( headl = .1 )
c
c draw vector
      call jsmove( x1, y1 )
      call jsline( x2, y2 )
c add arrow head - compute length using internal dimensions
      dx = xfac * ( x2 - x1 )
      dy = yfac * ( y2 - y1 )
      x = x2 + ( headl * ( dy - dx ) ) / xfac
      y = y2 - ( headl * ( dx + dy ) ) / yfac
      call jsline( x, y )
      x = x2 - ( headl * ( dx + dy ) ) / xfac
      y = y2 + ( headl * ( dx - dy ) ) / yfac
      call jsmove( x2, y2 )
      call jsline( x, y )
      return
      end
