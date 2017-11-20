      real function linint( x, y, n, xx )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c Linear interpolation scheme - works for unequally spaced abscissae
c   but they must be distinct and ordered
c
c calling arguments
      integer n
      real xx, x( n ), y( n )
c
c local variables
      integer i, nl, nu
      real p
c
c check that requested value is within tabulated range
      p = ( xx - x( 1 ) ) / ( x( n ) - x( 1 ) )
      if( ( p .lt. -1.e-4 ) .or. ( ( p - 1. ) .gt. 1.e-4 ) )then
c error message
        print *, 'range of x is', x( 1 ), ' to', x( n )
        print *, 'abscissa is', xx
        print *, 'number of ordinates', n
        call crash( 'LININT', 'Requested value outside range' )
      end if
c find straddling abscissae by bi-section - assumed in ascending order
      nl = 1
      nu = n
      do while ( nl + 1 .lt. nu )
        i = ( nl + nu ) / 2
        if( x( i ) .lt. xx )then
          nl = i
        else
          nu = i
        end if
      end do
      if( nl .ge. nu )call crash( 'LININT', 'Abscissae not ranked' )
c linear interpolation
      p = ( xx - x( nl ) ) / ( x( nu ) - x( nl ) )
      linint = ( 1. - p ) * y( nl ) + p * y( nu )
      return
      end

      real*8 function linint2( x, y, n, xx )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c Linear interpolation scheme - works for unequally spaced abscissae
c   but they must be distinct and ordered
c
c calling arguments
      integer n
      real*8 xx, x( n ), y( n )
c
c local variables
      integer i, nl, nu
      real*8 p
c
c check that requested value is within tabulated range
      p = ( xx - x( 1 ) ) / ( x( n ) - x( 1 ) )
      if( ( p .lt. -1.d-6 ) .or. ( ( p - 1. ) .gt. 1.d-6 ) )then
c error message
        print *, 'range of x is', x( 1 ), ' to', x( n )
        print *, 'abscissa is', xx
        print *, 'number of ordinates', n
        call crash( 'LININT2', 'Requested value outside range' )
      end if
c find straddling abscissae by bi-section - assumed in ascending order
      nl = 1
      nu = n
      do while ( nl + 1 .lt. nu )
        i = ( nl + nu ) / 2
        if( x( i ) .lt. xx )then
          nl = i
        else
          nu = i
        end if
      end do
      if( nl .ge. nu )call crash( 'LININT2', 'Abscissae not ranked' )
c linear interpolation
      p = ( xx - x( nl ) ) / ( x( nu ) - x( nl ) )
      linint2 = ( 1. - p ) * y( nl ) + p * y( nu )
      return
      end
