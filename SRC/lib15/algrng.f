      real function algrng( x, y, n, xx )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Lagrange interpolation scheme - equally spaced abscissae
c
c calling argumemts
      integer n
      real xx, x( n ), y( n )
c
c local variables
      integer i
      real p, pp
c
c check that requested value is within tabulated range
      p = ( xx - x( 1 ) ) / ( x( n ) - x( 1 ) )
      if( ( p .lt. -1.e-4 ) .or. ( ( p - 1. ) .gt. 1.e-4 ) )then
c error message
        print *, 'Requested value outside range - error in ALGRNG'
        print *, 'range of x is', x( 1 ), ' to', x( n )
        print *, 'abscissa is', xx
        print *, 'number of ordinates', n
        stop
      end if
c
      if( n .gt. 4 )then
        p = p * real( n - 1 )
        i = int( p )
        i = max( i, 2 )
        i = min( i, n - 3 )
        p = p - real( i )
c formula 25.2.15 from Abramowitz and Stegun
        pp = p * p
        algrng = p * ( pp - 1. ) * ( ( p - 2. ) * y( i - 1 )
     +                             + ( p + 2. ) * y( i + 3 ) ) / 24.
     +         - p * ( pp - 4. ) * ( ( p - 1. ) * y( i )
     +                             + ( p + 1. ) * y( i + 2 ) ) / 6.
     +         + .25 * ( pp - 1. ) * ( pp - 4. ) * y( i + 1 )
      else
c formula 25.2.13 from Abramowitz and Stegun
        if( n .eq. 4 )then
          p = 3. * p - 1.
          algrng = p * ( p - 1. ) * ( ( p + 1. ) * y( 4 )
     +                              - ( p - 2. ) * y( 1 ) ) / 6.
     +            + .5 * ( p + 1. ) * ( p - 2. ) * ( ( p - 1. ) * y( 2 )
     +                                                    - p * y( 3 ) )
c formula 25.2.11 from Abramowitz and Stegun
        else if( n .eq. 3 )then
          p = 2. * p - 1.
          algrng = .5 * p * ( ( p - 1. ) * y( 1 )
     +                      + ( 1. + p ) * y( 3 ) )
     +                      + ( 1. - p * p ) * y( 2 )
c linear interpolation
        else if( n .eq. 2 )then
          algrng = ( 1. - p ) * y( 1 ) + p * y( 2 )
        else
          print *, 'Impossible number of ordinates - error in ALGRNG'
          print *, 'number of ordinates', n
          stop
        end if
      end if
      return
      end

      real*8 function algrng2( x, y, n, xx )
c Lagrange interpolation scheme - equally spaced abscissae
c
c calling argumemts
      integer n
      real*8 xx, x( n ), y( n )
c
c local variables
      integer i
      real*8 p, pp
c
c check that requested value is within tabulated range
      p = ( xx - x( 1 ) ) / ( x( n ) - x( 1 ) )
      if( ( p .lt. -1.d-6 ) .or. ( ( p - 1.d0 ) .gt. 1.d-6 ) )then
c error message
        print *, 'Requested value outside range - error in ALGRNG2'
        print *, 'range of x is', x( 1 ), ' to', x( n )
        print *, 'abscissa is', xx
        print *, 'number of ordinates', n
        stop
      end if
c
      if( n .gt. 4 )then
        p = p * dble( n - 1 )
        i = int( p )
        i = max( i, 2 )
        i = min( i, n - 3 )
        p = p - dble( i )
c formula 25.2.15 from Abramowitz and Stegun
        pp = p * p
        algrng2 = p * ( pp - 1.d0 ) * ( ( p - 2.d0 ) * y( i - 1 )
     +                             + ( p + 2.d0 ) * y( i + 3 ) ) / 24.d0
     +          - p * ( pp - 4.d0 ) * ( ( p - 1.d0 ) * y( i )
     +                              + ( p + 1.d0 ) * y( i + 2 ) ) / 6.d0
     +          + .25d0 * ( pp - 1.d0 ) * ( pp - 4.d0 ) * y( i + 1 )
      else
c formula 25.2.13 from Abramowitz and Stegun
        if( n .eq. 4 )then
          p = 3.d0 * p - 1.d0
          algrng2 = p * ( p - 1.d0 ) * ( ( p + 1.d0 ) * y( 4 )
     +                               - ( p - 2.d0 ) * y( 1 ) ) / 6.d0
     +    + .5d0 * ( p + 1.d0 ) * ( p - 2.d0 ) * ( ( p - 1.d0 ) * y( 2 )
     +                                                    - p * y( 3 ) )
c formula 25.2.11 from Abramowitz and Stegun
        else if( n .eq. 3 )then
          p = 2.d0 * p - 1.d0
          algrng2 = .5d0 * p * ( ( p - 1.d0 ) * y( 1 )
     +                       + ( 1.d0 + p ) * y( 3 ) )
     +                       + ( 1.d0 - p * p ) * y( 2 )
c linear interpolation
        else if( n .eq. 2 )then
          algrng2 = ( 1.d0 - p ) * y( 1 ) + p * y( 2 )
        else
          print *, 'Impossible number of ordinates - error in ALGRNG2'
          print *, 'number of ordinates', n
          stop
        end if
      end if
      return
      end
