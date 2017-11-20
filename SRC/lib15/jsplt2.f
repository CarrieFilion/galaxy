      subroutine jsplt2( funct )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to draw the real*8 function in the specified window.
c   it is adaptive, and will space values more closely where the
c   slope changes rapidly
c
c calling argument
      real*8 funct
c
c common block
c
      include 'inc/jscmmn.f'
c
c local arrays
      real x( 13 ), y( 13 )
c
c local variables
      integer i, j, m, n
      parameter ( n = 50 )
      real dm, slopen, slopeo, xs, yp
c
      xs = ( sx2 - sx1 ) / real( 12 * n )
c calculate initial slope
      x( 1 ) = sx1
      y( 1 ) = funct( dble( x( 1 ) ) )
      yp = max( y( 1 ), sy1 )
      yp = min( yp, sy2 )
      call jsmove( x( 1 ), yp )
      x( 2 ) = x( 1 ) + xs
      y( 2 ) = funct( dble( x( 2 ) ) )
      slopeo = ( y( 2 ) - y( 1 ) ) / ( x( 2 ) - x( 1 ) )
c work over segments of curve
      do i = 1, n
        x( 13 ) = x( 1 ) + 12. * xs
        y( 13 ) = funct( dble( x( 13 ) ) )
c subdivide when slope changes substantially
        slopen = ( y( 13 ) - y( 1 ) ) / ( x( 13 ) - x( 1 ) )
        dm = abs( slopen - slopeo )
        m = 12
        if( dm .gt. .166666666 )m = int( 1. / dm ) + 1
        if( m .eq. 5 )m = 6
        do j = 1, 12, m
          if( j .ne. 1 )then
            x( j ) = x( 1 ) + real( j - 1 ) * xs
            y( j ) = funct( dble( x( j ) ) )
            slopen = ( y( 13 ) - y( j ) ) / ( x( 13 ) - x( j ) )
          end if
        end do
c plot calculated points
        do j = 1, 13, m
          if( j .ne. 1 )call jsline( x( j ), y( j ) )
        end do
c set values for next segment
        slopeo = slopen
        x( 1 ) = x( 13 )
        y( 1 ) = y( 13 )
      end do
      return
      end
