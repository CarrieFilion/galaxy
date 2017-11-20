      subroutine intspl( m, x, y, k, c )
c determines the least-squares solution for cubic-spline interpolation
c   giving equal weight to each data point
c follows the algorithm in the NAG source code
      implicit none
c
c calling arguments
      integer m
      real c( m + 4 ), k( m + 4 ), x( m ), y( m )
c
c local workspace
      real, allocatable :: w1( : ), w2( :, : )
c
c local variables
      integer i, ij, iu, j, j0, jl, jrev, l, l4, l1, lu, n7, np
      integer np3, nm1, r
      real a1, a2, b1, b2, b3, b4, d, d4, d5, d6, d7, d8, d9
      real dp, e2, e3, e4, e5, n1, n2, n3, rt, s
c
c check that the value of m is reasonable
      if( m .lt. 4 )call crash( 'INTSPL', 'Too few abscissae' )
c set knots for an interpolating spline
      do i = 4, m
        k( i ) = x( i - 2 )
      end do
      do i = 1, 4
        k( i ) = x( 1 )
        k( m + i ) = x( m )
      end do
c allocate work arrays
      n7 = m + 4
      allocate ( w1( m ) )
      allocate ( w2( 4, n7 ) )
c check that the abscissae are ordered and copy only distinct abscissae
c   to the array w1
      w1( 1 ) = x( 1 )
      j = 2
      do i = 2, m
        if( x( i ) .lt.
     +      w1( j - 1 ) )call crash( 'INTSPL', 'Abscissae not ordered' )
        if( x( i ) .ne. w1( j - 1 ) )then
          w1( j ) = x( i )
          j = j + 1
        end if
      end do
      r = j - 1
c check that there are sufficient distinct abscissae for the number of knots
      np = n7 - 7
      np3 = np + 3
      if(
     +  r .lt. np3 )call crash( 'INTSPL', 'Too few distinct abscissae' )
c check the first  s  and the last  s  Schoenberg-Whitney
c   conditions ( s = min( np - 1, 4 ) )
      do j = 1, 4
        if( j .ge. np )go to 2
        i = np3 - j + 1
        l = r - j + 1
        if( ( w1( j ) .ge. k( j + 4 ) ) .or. ( k( i ) .ge. w1( l ) )
     +               )call crash( 'INTSPL', 'Matrix probably singular' )
      end do
c check all the remaining Schoenberg-Whitney conditions
      nm1 = np - 1
      if( np .gt. 5 )then
        r = r - 4
        i = 4
        do j = 5, nm1
          i = i + 1
          do while ( w1( i ) .le. k( j ) )
            i = i + 1
          end do
          if( i .gt. r .or. w1( i ) .ge. k( j + 4 )
     +               )call crash( 'INTSPL', 'Matrix probably singular' )
        end do
      end if
c initialize the tri-diagonal matrix
    2 do i = 1, np3
        do l = 1, 4
          w2( l, i ) = 0.0
        end do
        c( i ) = 0.0
      end do
      j = 0
      j0 = 0
      do i = 1, m
c find the knots that bracket each x( i ) in turn
        do while ( x( i ) .ge. k( j + 4 ) .and. j .le. nm1 )
          j = j + 1
        end do
        if( j .ne. j0 )then
c set constants for this j
          d4 = 1. / ( k( j + 4 ) - k( j + 1 ) )
          d5 = 1. / ( k( j + 5 ) - k( j + 2 ) )
          d6 = 1. / ( k( j + 6 ) - k( j + 3 ) )
          d7 = 1. / ( k( j + 4 ) - k( j + 2 ) )
          d8 = 1. / ( k( j + 5 ) - k( j + 3 ) )
          d9 = 1. / ( k( j + 4 ) - k( j + 3 ) )
          j0 = j
        end if
c compute and store in  w1(l ) (l = 1, 2, 3, 4 )  the values of the
c   four normalized cubic B-splines which are non-zero at x = x( i )
        e5 = k( j + 5 ) - x( i )
        e4 = k( j + 4 ) - x( i )
        e3 = x( i ) - k( j + 3 )
        e2 = x( i ) - k( j + 2 )
        n1 = d9
        n2 = e3 * n1 * d8
        n1 = e4 * n1 * d7
        n3 = e3 * n2 * d6
        n2 = ( e2 * n1 + e5 * n2 ) * d5
        n1 = e4 * n1 * d4
        w1( 4 ) = e3 * n3
        w1( 3 ) = e2 * n2 + ( k( j + 6 ) - x( i ) ) * n3
        w1( 2 ) = ( x( i ) - k( j + 1 ) ) * n1 + e5 * n2
        w1( 1 ) = e4 * n1
        b2 = y( i )
c rotate this row into the band triangular system using plane rotations
        do l1 = 1, 4
          l = l1 - 1
          rt = w1( l1 )
          if( rt .ne. 0. )then
            jl = j + l
            l4 = 4 - l
            d = w2( 1, jl )
            if( abs( rt ) .ge. d )then
              dp = abs( rt ) * sqrt( 1. + ( d / rt )**2 )
            else
              dp = d * sqrt( 1. + ( rt / d )**2 )
            end if
            w2( 1, jl ) = dp
            b3 = d / dp
            b4 = rt / dp
            if( l4 .ge. 2 )then
              do iu = 2, l4
                lu = l + iu
                a1 = w2( iu, jl )
                a2 = w1( lu )
                w2( iu, jl ) = b3 * a1 + b4 * a2
                w1( lu ) = b3 * a2 - b4 * a1
              end do
            end if
            b1 = c( jl )
            c( jl ) = b3 * b1 + b4 * b2
            b2 = b3 * b2 - b4 * b1
          end if
        end do
      end do
c solve the tri-diagonal system for the B-spline coefficients
      l = -1
      do jrev = 1, np3
        j = np3 - jrev + 1
        d = w2( 1, j )
        if( d .eq. 0. )call crash( 'INTSPL', 'Matrix is singular' )
        if( l .lt. 3 )l = l + 1
        s = c( j )
        if( l .gt. 0 )then
          do i = 1, l
            ij = i + j
            s = s - w2( i + 1, j ) * c( ij )
          end do
        end if
        c( j ) = s / d
      end do
      return
      end

      subroutine intspl2( m, x, y, k, c )
c determines the least-squares solution for cubic-spline interpolation
c   giving equal weight to each data point
      implicit none
c
c calling arguments
      integer m
      real*8 c( m + 4 ), k( m + 4 ), x( m ), y( m )
c
c local workspace
      real*8, allocatable :: w1( : ), w2( :, : )
c
c local variables
      integer i, ij, iu, j, j0, jl, jrev, l, l4, l1, lu, n7, np
      integer np3, nm1, r
      real*8 a1, a2, b1, b2, b3, b4, d, d4, d5, d6, d7, d8, d9
      real*8 dp, e2, e3, e4, e5, n1, n2, n3, rt, s
c
c check that the value of m is reasonable
      if( m .lt. 4 )call crash( 'INTSPL2', 'Too few abscissae' )
c set knots for an interpolating spline
      do i = 4, m
        k( i ) = x( i - 2 )
      end do
      do i = 1, 4
        k( i ) = x( 1 )
        k( m + i ) = x( m )
      end do
c allocate work arrays
      n7 = m + 4
      allocate ( w1( m ) )
      allocate ( w2( 4, n7 ) )
c check that the abscissae are ordered and copy only distinct abscissae
c   to the array w1
      w1( 1 ) = x( 1 )
      j = 2
      do i = 2, m
        if( x( i ) .lt.
     +     w1( j - 1 ) )call crash( 'INTSPL2', 'Abscissae not ordered' )
        if( x( i ) .ne. w1( j - 1 ) )then
          w1( j ) = x( i )
          j = j + 1
        end if
      end do
      r = j - 1
c check that there are sufficient distinct abscissae for the number of knots
      np = n7 - 7
      np3 = np + 3
      if(
     + r .lt. np3 )call crash( 'INTSPL2', 'Too few distinct abscissae' )
c check the first  s  and the last  s  Schoenberg-Whitney
c   conditions ( s = min( np - 1, 4 ) )
      do j = 1, 4
        if( j .ge. np )go to 2
        i = np3 - j + 1
        l = r - j + 1
        if( ( w1( j ) .ge. k( j + 4 ) ) .or. ( k( i ) .ge. w1( l ) )
     +           )call crash( 'INTSPL2', 'Matrix is probably singular' )
      end do
c check all the remaining Schoenberg-Whitney conditions
      nm1 = np - 1
      if( np .gt. 5 )then
        r = r - 4
        i = 4
        do j = 5, nm1
          i = i + 1
          do while ( w1( i ) .le. k( j ) )
            i = i + 1
          end do
          if( i .gt. r .or. w1( i ) .ge. k( j + 4 )
     +           )call crash( 'INTSPL2', 'Matrix is probably singular' )
        end do
      end if
c initialize the tri-diagonal matrix
    2 do i = 1, np3
        do l = 1, 4
          w2( l, i ) = 0.0d0
        end do
        c( i ) = 0.0d0
      end do
      j = 0
      j0 = 0
      do i = 1, m
c find the knots that bracket each x( i ) in turn
        do while ( x( i ) .ge. k( j + 4 ) .and. j .le. nm1 )
          j = j + 1
        end do
        if( j .ne. j0 )then
c set constants for this j
          d4 = 1.d0 / ( k( j + 4 ) - k( j + 1 ) )
          d5 = 1.d0 / ( k( j + 5 ) - k( j + 2 ) )
          d6 = 1.d0 / ( k( j + 6 ) - k( j + 3 ) )
          d7 = 1.d0 / ( k( j + 4 ) - k( j + 2 ) )
          d8 = 1.d0 / ( k( j + 5 ) - k( j + 3 ) )
          d9 = 1.d0 / ( k( j + 4 ) - k( j + 3 ) )
          j0 = j
        end if
c compute and store in  w1(l ) (l = 1, 2, 3, 4 )  the values of the
c   four normalized cubic B-splines which are non-zero at x = x( i )
        e5 = k( j + 5 ) - x( i )
        e4 = k( j + 4 ) - x( i )
        e3 = x( i ) - k( j + 3 )
        e2 = x( i ) - k( j + 2 )
        n1 = d9
        n2 = e3 * n1 * d8
        n1 = e4 * n1 * d7
        n3 = e3 * n2 * d6
        n2 = ( e2 * n1 + e5 * n2 ) * d5
        n1 = e4 * n1 * d4
        w1( 4 ) = e3 * n3
        w1( 3 ) = e2 * n2 + ( k( j + 6 ) - x( i ) ) * n3
        w1( 2 ) = ( x( i ) - k( j + 1 ) ) * n1 + e5 * n2
        w1( 1 ) = e4 * n1
        b2 = y( i )
c rotate this row into the band triangular system using plane rotations
        do l1 = 1, 4
          l = l1 - 1
          rt = w1( l1 )
          if( rt .ne. 0.d0 )then
            jl = j + l
            l4 = 4 - l
            d = w2( 1, jl )
            if( abs( rt ) .ge. d )then
              dp = abs( rt ) * sqrt( 1.d0 + ( d / rt )**2 )
            else
              dp = d * sqrt( 1.d0 + ( rt / d )**2 )
            end if
            w2( 1, jl ) = dp
            b3 = d / dp
            b4 = rt / dp
            if( l4 .ge. 2 )then
              do iu = 2, l4
                lu = l + iu
                a1 = w2( iu, jl )
                a2 = w1( lu )
                w2( iu, jl ) = b3 * a1 + b4 * a2
                w1( lu ) = b3 * a2 - b4 * a1
              end do
            end if
            b1 = c( jl )
            c( jl ) = b3 * b1 + b4 * b2
            b2 = b3 * b2 - b4 * b1
          end if
        end do
      end do
c solve the tri-diagonal system for the B-spline coefficients
      l = -1
      do jrev = 1, np3
        j = np3 - jrev + 1
        d = w2( 1, j )
        if( d .eq. 0.d0 )call crash( 'INTSPL2', 'Matrix is singular' )
        if( l .lt. 3 )l = l + 1
        s = c( j )
        if( l .gt. 0 )then
          do i = 1, l
            ij = i + j
            s = s - w2( i + 1, j ) * c( ij )
          end do
        end if
        c( j ) = s / d
      end do
      return
      end
