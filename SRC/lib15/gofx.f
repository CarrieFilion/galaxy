      real*8 function gofx( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c evaluates Kalnajs' expression (27) on p 755 of Ap J 205 (1976)
c  needed to evaluate the distribution function of type e**(m-1)*g(x)
c
c calling argument
      real*8 x
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      common / kalndf / xx, m
      integer m
      real*8 xx
c
      integer itab, iuseg
      parameter ( itab = 501 )
      common / gtab / xmax, xt( itab ), gt( itab ), iuseg
      real*8 xmax, xt, gt
c
      include 'inc/model.f'
c
c externals
      real*8 algrng2, deriv2, gofxkt, quad_Pat, phitot, vcirc, ymtau
      external gfkaln, ymtau
c
c local variables
      integer ifail, ix
      real*8 aa, bb, dfm, epsr
      include 'inc/pi.f'
c
      dfm = dfcns( 1, icmp )
      m = dfm + .5
      if( iuseg .ne. m )then
c compute xmax
        aa = 1.01 * rmax
        bb = phitot( aa ) + .5 * vcirc( aa )**2
        xmax = -aa * vcirc( aa ) * sqrt( -2. * bb ) / ( phi0 * rstar )
c construct table of values
        do ix = 1, itab
          xx = xmax * dble( ix - 1 ) / dble( itab - 1 )
          xt( ix ) = xx
          if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'KT  ' ) )then
c analytic expression
            gofx = gofxkt( xx )
          else
c special value of argument
            if( xx .eq. 0. )then
              gofx = 0
              if( ctype( 1 ) .eq. 'SOFT' )gofx =
     +                                      dfm / ( -2. * pi**2 * phi0 )
              if( ctype( 1 ) .eq. 'ISOC' )gofx = dfm / ( 3. * pi**2 )
              if( gofx .eq. 0. )then
                print *, 'Unrecognised model - error in GOFX'
                stop
              end if
            else
c general expression for arbitrary disc - eq. (27) of Kalnajs
              aa = 0
              bb = 1.
              epsr = 1.e-6
              ifail = 0
              gofx = quad_Pat( gfkaln, aa, bb, epsr, ifail )
              gofx = gofx + xx * deriv2( xx, ymtau )
     +                    - .5 * dfm * ( dfm - 1. ) * ymtau( xx )
c multiplicative factor
              gofx = ( 2. / xx )**m * gofx / ( -2. * pi * phi0 )
            end if
          end if
          if( gofx .lt. 0. )then
            print *, ix, xx, gofx
            call crash( 'GOFX', 'gone negative' )
c            gofx = 10.**gt( ix - 1 )
          end if
          gt( ix ) = log10( gofx )
        end do
        iuseg = m
      end if
c Lagrange interpolation
      if( abs( x ) .gt. xmax )then
        print *, 'x = ', x
        call crash( 'GOFX', 'Impossible argument' )
      end if
      gofx = algrng2( xt, gt, itab, abs( x ) )
      gofx = 10.**gofx
      return
      end
