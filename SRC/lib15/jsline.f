      subroutine jsline( x, y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to draw a line from the current pen position to the point (x,y)
c   using the pre-specified line style and weight
c
c calling arguments
      real x, y
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local variables
      integer i
      real cth, d, dist, sth, vec, xn, yn
c
      if( dash )then
c draw dashed line - compute length and direction
        d = sqrt( ( ( x - xlast ) * xfac )**2 +
     +            ( ( y - ylast ) * yfac )**2 )
        if( d .gt. 0. )then
          cth = ( x - xlast ) * xfac / d
          sth = ( y - ylast ) * yfac / d
c initialize
          xn = xlast
          yn = ylast
          dist = d
          if( partd .ne. 0. )then
c remaining length to draw
            vec = 0
            do i = 1, iseg
              vec = vec + segl( i )
            end do
            vec = vec - partd
          else
c draw a full segment
            vec = segl( iseg )
          end if
c draw line segments
          do while ( dist .gt. 0. )
            vec = min( vec, dist )
            xn = xn + vec * cth / xfac
            yn = yn + vec * sth / yfac
            if( mod( iseg, 2 ) .eq. 0 )then
              call jspenu( xn, yn )
            else
              call jspend( xn, yn )
            end if
c deal with roundoff issue
            if( vec .eq. dist )then
              dist = 0
            else
              dist = dist - vec
              iseg = mod( iseg, 4 ) + 1
              vec = segl( iseg )
            end if
          end do
c store new partd
          partd = mod( partd + d, total )
        end if
      else
c draw solid line
        call jspend( x, y )
      end if
c save new position
      xlast = x
      ylast = y
      return
      end
