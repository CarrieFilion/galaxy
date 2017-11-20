      real*8 function poveda( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c density as a function of r in a spherical r^(1/4) model by the Abel
c   integral form derived by Poveda - ref: Young (1976, AJ 81, 807)
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
      common / povedal / b, aj
      real*8 aj, b
c
c externals
      external povfun
      real*8 algrng2, quad_inf
c
c local arrays
      integer ntab
      parameter ( ntab = 50 )
      real*8 x( ntab ), y( ntab )
      save x, y
c
c local variables
      integer ier, inf, iuse, j
      real*8 bound, epsa, epsr, s, tab, young
      save iuse
      include 'inc/pi.f'
c
      data iuse / 0 /
c
      poveda = 0
      if( r .gt. 0. )then
        if( iuse .eq. 0 )then
c determine b - Young's equation (4) with alpha = 1 and f = 0.5
          b = 7
          s = 8
          inf = 1
          epsr = 1.e-20
          ier = 1
c find zero of function
          do while ( inf .ne. 0 )
            call fndzro( b, s, bound, epsr, 0, x, inf, ier )
            if( inf .ne. 0 )then
              tab = exp( -b )
              bound = 0.5 - tab
              do j = 1, 7
                tab = tab * b / dble( j )
                bound = bound - tab
              end do
            end if
          end do
          if( master )print *, 'Building a table of bulge densities'
          do j = 1, ntab
c equal spacing in log(r) over range: exp(-18.5) < r/R_e < exp(6)=403.4
            s = j + 12 - ntab
            s = .5 * s
            x( j ) = s
            aj = exp( .25 * s )
c Poveda's numerical integral from 1 to infinity - Young's eq (7)
            bound = 1
            inf = 1
            epsa = 1.e-20
            epsr = 1.e-8
            tab = quad_inf( povfun, bound, inf, epsa, epsr, ier )
            if( ier .gt. 0 )then
              print *, 'ier =', ier, ' from QUAD_INF'
              call crash( 'POVEDA', 'QUADPACK error' )
            end if
c store deviation from Young's approximate eq (33) - also without 1/(2j**3)
            young = sqrt( pi / ( 8. * b * aj ) ) * exp( -b * aj )
            y( j ) = tab / young
          end do
          iuse = 1
          if( master )print *, 'Table ready'
        end if
c table in equal steps in log( r )
        aj = r**.25
        s = log( r )
        if( s .lt. x( 1 ) )then
c Young's approximate eq (29) for very small r
          poveda = 0.2409903 * exp( -b * aj ) / aj**3
        else
c look up tabulated value - good to r = 400 * R_e and beyond
          s = min( s, x( ntab ) )
          tab = algrng2( x, y, ntab, s )
c factor in approximate expression
          young = sqrt( pi / ( 8. * b * aj ) ) * exp( -b * aj )
     +                                                  / ( 2. * aj**3 )
          poveda = tab * young
        end if
c normalise to unit total mass - the constant is 2 b**9 / ( 8! pi**2 )
        poveda = 461.302165 * poveda
      end if
      return
      end
