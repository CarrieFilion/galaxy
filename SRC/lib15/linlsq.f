      subroutine linlsq( xval, yval, n, slope, bint, slerr, berr )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c linear regression to find least-squares fit
c
c calling arguments
      integer n
      real berr, bint, slerr, slope, xval( n ), yval( n )
c
c local variables
      integer i
      real delta, t, x, xx, xy, y
c
c find slope and intercept
      x = 0.
      xx = 0.
      y = 0.
      xy = 0.
      do i = 1, n
        t = xval( i )
        x = x + t
        xx = xx + t * t
        xy = xy + t * yval( i )
        y = y + yval( i )
      end do
      t = n
      slope = ( t * xy - x * y ) / ( t * xx - x * x )
      bint = ( y * xx - xy * x ) / ( t * xx - x * x )
c determine precision of slope
      slerr = 0.
      do i = 1, n
        delta = slope * xval( i ) + bint - yval( i )
        slerr = slerr + delta * delta
      end do
      slerr = slerr / ( t * xx - x * x )
      slerr = abs( slerr )
      slerr = sqrt( slerr )
c code for intercept error not written
      return
      end
