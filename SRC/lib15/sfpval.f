      subroutine sfpval( r )
      use aarrays
c  Copyright (C) 2014, Jerry Sellwood
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
      integer i, j
      real*8 p, pm1, pm2, ppm1, ppm4, pp1, pp2, w1, w2, w3
      parameter ( w1 = 1.d0 / 24.d0, w2 = 1.d0 / 6.d0, w3 = .25d0 )
c
c tabulated function values
      p = -1
      if( basset .eq. 'ablj' )p = r / maxr
      if( basset .eq. 'bess' )p = r / sfprad( lsfpr )
c check that requested value is within tabulated range
      if( ( p .lt. 0.d0 ) .or. ( p .gt. 1.d0 ) )then
        print *, 'abscissa is', r
        call crash( 'SFPVAL', 'Requested value outside range' )
      end if
c
      p = p * dble( lsfpr - 1 )
      i = p
      i = max( i, 2 )
      i = min( i, lsfpr - 3 )
      p = p - dble( i )
c formula 25.2.15 from Abramowitz & Stegun
      pm1 = p - 1.
      pm2 = p - 2.
      pp1 = p + 1.
      pp2 = p + 2.
      ppm1 = p * p - 1.
      ppm4 = p * p - 4.
      do j = 1, lastf
        funvals( j ) = p * ppm1 * ( pm2 * sfplst( 1, i - 1 )
     +                            + pp2 * sfplst( 1, i + 3 ) ) * w1
     +               - p * ppm4 * ( pm1 * sfplst( 1, i )
     +                            + pp1 * sfplst( 1, i + 2 ) ) * w2
     +                    + ppm1 * ppm4 * sfplst( 1, i + 1 )   * w3
        i = i + lsfpr
      end do
      return
      end
