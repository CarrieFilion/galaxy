      real function splval( x, n, k, c )
c evaluates a cubic spline from its B-spline representation
c   follows the algorithm in the NAG source code
      implicit none
c
c calling arguments
      integer n
      real c( n ), k( n ), x
c
c local variables
      integer i, j, l
      real a1, a2, a3, d2, d3, d4, d5
c
c check calling arguments
      if( n .lt. 8 )call crash( 'SPLVAL', 'Too few knots' )
      if( x .lt. k( 4 ) .or. x .gt. k( n - 3 )
     +         )call crash( 'SPLVAL', 'Argument outside allowed range' )
c find pair of bracketing knots
      i = 0
      j = n - 7
      do while ( j - i .gt. 1 )
        l = ( i + j ) / 2
        if( x .lt. k( l + 4 ) )then
          j = l
        else
          i = l
        end if
      end do
c de Boor's method of convex combinations
      d2 = x - k( j + 2 )
      d3 = x - k( j + 3 )
      d4 = k( j + 4 ) - x
      d5 = k( j + 5 ) - x
      a1 = ( ( x - k( j + 1 ) ) * c( j + 1 ) + d4 * c( j ) ) /
     +                                       ( k( j + 4 ) - k( j + 1 ) )
      a2 = ( d2 * c( j + 2 ) + d5 * c( j + 1 ) ) /
     +                                       ( k( j + 5 ) - k( j + 2 ) )
      a3 = ( d3 * c( j + 3 ) + ( k( j + 6 ) - x ) * c( j + 2 ) ) /
     +                                       ( k( j + 6 ) - k( j + 3 ) )
      a1 = ( d2 * a2 + d4 * a1 ) / ( k( j + 4 ) - k( j + 2 ) )
      a2 = ( d3 * a3 + d5 * a2 ) / ( k( j + 5 ) - k( j + 3 ) )
c set function value
      splval = ( d3 * a2 + d4 * a1 ) / ( k( j + 4 ) - k( j + 3 ) )
      return
      end

      real*8 function splval2( x, n, k, c )
c evaluates a cubic spline from its B-spline representation
      implicit none
c
c calling arguments
      integer n
      real*8 c( n ), k( n ), x
c
c local variables
      integer i, j, l
      real*8 a1, a2, a3, d2, d3, d4, d5
c
c check calling arguments
      if( n .lt. 8 )call crash( 'SPLVAL2', 'Too few knots' )
      if( x .lt. k( 4 ) .or. x .gt. k( n - 3 )
     +        )call crash( 'SPLVAL2', 'Argument outside allowed range' )
c find pair of bracketing knots
      i = 0
      j = n - 7
      do while ( j - i .gt. 1 )
        l = ( i + j ) / 2
        if( x .lt. k( l + 4 ) )then
          j = l
        else
          i = l
        end if
      end do
c de Boor's method of convex combinations
      d2 = x - k( j + 2 )
      d3 = x - k( j + 3 )
      d4 = k( j + 4 ) - x
      d5 = k( j + 5 ) - x
      a1 = ( ( x - k( j + 1 ) ) * c( j + 1 ) + d4 * c( j ) ) /
     +                                       ( k( j + 4 ) - k( j + 1 ) )
      a2 = ( d2 * c( j + 2 ) + d5 * c( j + 1 ) ) /
     +                                       ( k( j + 5 ) - k( j + 2 ) )
      a3 = ( d3 * c( j + 3 ) + ( k( j + 6 ) - x ) * c( j + 2 ) ) /
     +                                       ( k( j + 6 ) - k( j + 3 ) )
      a1 = ( d2 * a2 + d4 * a1 ) / ( k( j + 4 ) - k( j + 2 ) )
      a2 = ( d3 * a3 + d5 * a2 ) / ( k( j + 5 ) - k( j + 3 ) )
c set function value
      splval2 = ( d3 * a2 + d4 * a1 ) / ( k( j + 4 ) - k( j + 3 ) )
      return
      end
