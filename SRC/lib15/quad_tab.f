      real*8 function quad_tabc( x, y, n, er )
      implicit none
c front end routine for numerical integration of a function specified by
c   a table of values only - it just checks the calling arguments and then
c   calls the external that does the work
c
c arguments:  n     the number of values, which must be greater than 4
c             x, y  the abscissae and ordinates - unchanged on exit
c             er    an estimate of the uncertainty
c
c calling arguments
      integer n
      real*8 er, x( n ), y( n )
c
c external
      real*8 quad_tab
c
c local variables
      integer i
      real*8 d2, d3
c
c check arguments
      if( n .lt. 4 )call crash( 'QUAD_TAB', 'Too few abscissae' )
      d2 = x( 2 ) - x( 1 )
      do i = 3, n
        d3 = x( i ) - x( i - 1 )
        if( d2 * d3 .lt. 0.d0 )then
          call crash( 'QUAD_TAB', 'Abscissae not ordered' )
        else if( d2 * d3 .eq. 0.d0 )then
          call crash( 'QUAD_TAB', 'Abscissae not distinct' )
        end if
        d2 = d3
      end do
      quad_tabc = quad_tab( x, y, n, er )
      return
      end

      real*8 function quad_tab( x, y, n, er )
      implicit none
c returns an estimate of the integral of a function specified by a
c   ordinates, y(i), at abscissae x(i).  The integral is over the
c   entire interval x(1) -> x(n).  The n abscissae need not be equally
c   spaced, but should be distinct and in ascending or descending order.
c The method is due to Gill and Miller and follows the algorithm given
c   in their paper (Computer Journal (1972) ,15 p80)
c
c arguments:  n     the number of values, which must be greater than 4
c             x, y  the abscissae and ordinates - unchanged on exit
c             er    an estimate of the uncertainty
c
c calling arguments
      integer n
      real*8 er, x( n ), y( n )
c
c local variables
      integer i
      real*8 a, c, d1, d2, d3, d4, m1, m2, m3, r1, r2, r3, r4, s
c
c first interval
      d2 = x( n ) - x( n - 1 )
      m3 = ( y( 2 ) - y( 1 ) ) / d2
      d3 = x( 3 ) - x( 2 )
      m1 = ( y( 3 ) - y( 2 ) ) / d3
      d1 = d2 + d3
      m2 = ( m1 - m3 ) / d1
      d4 = x( 4 ) - x( 3 )
      r1 = ( y( 4 ) - y( 3 ) ) / d4
      r2 = ( r1 - m1 ) / ( d4 + d3 )
      d1 = d1 + d4
      r3 = ( r2 - m2 ) / d1
      a = d2 * ( y( 1 ) + d2 * ( .5d0 * m3 -
     +          d2 * ( m2 / 6.d0 - ( d2 + 2.d0 * d3 ) * r3 / 12.d0 ) ) )
      s = -d2**3 *
     +      ( d2 * ( 3.d0 * d2 + 5.d0 * d4 ) + 10.d0 * d3 * d1 ) / 60.d0
      r4 = 0.d0
      er = 0.d0
c main range
      do i = 3, n - 1
        a = a + d3 * ( .5d0 * ( y( i ) + y( i - 1 ) ) -
     +                d3 * d3 * ( m2 + r2 + ( d2 - d4 ) * r3 ) / 12.d0 )
        c = d3**3 * ( 2.d0 * d3 * d3 +
     +          5.d0 * ( d3 * ( d4 + d2 ) + 2.d0 * d4 * d2 ) ) / 120.d0
        er = er + ( c + s ) * r4
        if( i .eq. 3 )then
          s = s + 2.d0 * c
        else
          s = c
        end if
        if( i .ne. n - 1 )then
          d1 = d2
          d2 = d3
          d3 = d4
          m1 = r1
          m2 = r2
          m3 = r3
          d4 = x( i + 2 ) - x( i + 1 )
          r1 = ( y( i + 2 ) - y( i + 1 ) ) / d4
          r4 = d4 + d3
          r2 = ( r1 - m1 ) / r4
          r4 = r4 + d2
          r3 = ( r2 - m2 ) / r4
          r4 = r4 + d1
          r4 = ( r3 - m3 ) / r4
        end if
      end do
c last interval
      a = a + d4 * ( y( n ) - d4 * ( .5d0 * r1 +
     +          d4 * ( r2 / 6.d0 + ( 2.d0 * d3 + d4 ) * r3 / 12.d0 ) ) )
c error estimate
      er = er - d4**3 * r4 * ( d4 * ( 3.d0 * d4 + 5.d0 * d2 ) +
     +                  10.d0 * d3 * ( d2 + d3 + d4 ) ) / 60.d0 + s * r4
c add error to be consistent with NAG
      quad_tab = a + er
      return
      end
