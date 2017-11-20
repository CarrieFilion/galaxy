      subroutine yrange( xmin, xmax, ymin, ymax, fun )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c routine to exploer the values of the given function "fun" over the
c   given range of abscissae in order to determine a suitable range
c   of ordinates for plotting
c
c calling arguments
      real xmin, xmax, ymin, ymax
      real*8 fun
c
c external
      real roundup
c
c local variables
      integer i
      real x, y
c
      ymin = 1.e+10
      ymax = -ymin
      do i = 1, 101
        x = .01 * real( i - 1 )
        x = ( 1. - x ) * xmin + x * xmax
        y = fun( dble( x ) )
        ymax = max( y, ymax )
        ymin = min( y, ymin )
      end do
      if( ( ymax .lt. 0. ) .and. ( -ymax .lt. -.5 * ymin ) )then
        ymax = 0
      else
        if( ( ymin .gt. 0. ) .and. ( ymin .lt. .3 * ymax ) )ymin = 0
      end if
      if( ymax .gt. 0. )ymax = roundup( ymax )
      if( ymin .gt. 0. )then
        ymin = roundup( ymin )
      else
        ymin = -roundup( -ymin )
      end if
      return
      end
