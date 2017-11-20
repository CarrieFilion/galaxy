      subroutine jscirc( x, y, r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to draw a "circle" radius r and center at (x,y)
c   the circle will in fact be an ellipse unless equal scaling is
c   used in the x- and y-directions
c
c calling arguments
      real r, x, y
c
c local variables
      integer i, n
      real pi, th, xx, yy
      parameter ( pi = 3.1415926535898, n = 100 )
c
      call jsmove( x + r, y )
      do i = 1, n
        th = 2. * pi * real( i ) / real( n )
        xx = r * cos( th ) + x
        yy = r * sin( th ) + y
        call jsline( xx, yy )
      end do
      return
      end
