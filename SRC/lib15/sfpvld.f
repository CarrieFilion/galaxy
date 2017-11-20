      subroutine sfpvld( r )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c Lagrange interpolation scheme - equally spaced abscissae
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
c local variables
      integer i, j, n
      real*8 p, pm1, pm2, ppm1, ppm4, pp1, pp2, w1, w2, w3
      parameter ( w1 = 1.d0 / 24.d0, w2 = 1.d0 / 6.d0, w3 = .25d0 )
c
c      if( ( basset .eq. 'lgsp' ) .and. ( r .lt. sfprad( 3 ) ) )then
      if( basset .eq. 'lgsp' )then
c ignore interpolation
        do j = 1, lastf
          n = nsel( j ) + 1
          p = real( n - maxn / 2 ) * deltak * log( r )
          funvals( j ) = cos( p ) / sqrt( r )
          dervals( j ) = sin( p ) / sqrt( r )
        end do
      else
        p = -1
        if( basset .eq. 'ablj' )p = r / maxr
c        if( ( basset .eq. 'bess' ) .or.
c     +      ( basset .eq. 'lgsp' ) )p = r / sfprad( lsfpr )
        if( basset .eq. 'bess' )p = r / sfprad( lsfpr )
c check that requested value is within 0 <= r <= 1
        if( ( p .lt. 0.d0 ) .or. ( p - 1.d0 .gt. 1.d-5 ) )then
          print *, 'Abscissa is', r
          call crash( 'SFPVLD', 'Requested value outside range' )
        end if
c
        p = p * dble( lsfpr - 1 )
        i = p
        i = max( i, 2 )
        i = min( i, lsfpr - 3 )
        p = p - dble( i )
c formula 25.2.15 from Abramowitz and Stegun
        pm1 = p - 1.
        pm2 = p - 2.
        pp1 = p + 1.
        pp2 = p + 2.
        ppm1 = p * p - 1.
        ppm4 = p * p - 4.
        do j = 1, lastf
          funvals( j ) = p * ppm1 * ( pm2 * sfplst( 1, i - 1 )
     +                              + pp2 * sfplst( 1, i + 3 ) ) * w1
     +                 - p * ppm4 * ( pm1 * sfplst( 1, i )
     +                              + pp1 * sfplst( 1, i + 2 ) ) * w2
     +                      + ppm1 * ppm4 * sfplst( 1, i + 1 )   * w3
          dervals( j ) = p * ppm1 * ( pm2 * sfplst( 2, i - 1 )
     +                              + pp2 * sfplst( 2, i + 3 ) ) * w1
     +                 - p * ppm4 * ( pm1 * sfplst( 2, i )
     +                              + pp1 * sfplst( 2, i + 2 ) ) * w2
     +                      + ppm1 * ppm4 * sfplst( 2, i + 1 )   * w3
          i = i + lsfpr
        end do
      end if
      return
      end
