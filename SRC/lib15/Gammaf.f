      real*8 function Gammaf( x, ifail )
c returns the value of the Gamma function for the input argument x
c   follows the algorithm in the NAG source code
c   the 2nd parameter, ifail, is redundant
c
c calling arguments
      integer ifail
      real*8 x
c
c local variables
      integer i, m
      real*8 g, gbig, t, xbig, xminv, xsmall, y
c
      data xsmall / 1.d-17 /
      data gbig, xbig, xminv / 4.3d+304, 170.d0, 2.23d-308 /
c
c     xbig = largest x such that  gamma(x) .lt. maxreal
c        and  1.0/gamma(x+1.0) .gt. minreal (rounded down to an integer)
c     gbig = gamma(xbig)
c     xminv = max(1.0/maxreal,minreal)  (rounded up)
c
      t = abs( x )
      m = t
c large argument
      if( t .gt. xbig )then
        if( x .gt. 0.d0 )then
          Gammaf = gbig
        else
          Gammaf = 0
        end if
        if( ifail .eq. 0 )call crash( 'GAMMAF', 'Argument too large' )
        ifail = 1
c small argument
      else if( t .le. xsmall )then
        t = max( t, xminv )
        Gammaf = sign( 1.d0 / t, x )
        if( ifail .eq. 0 )call crash( 'GAMMAF', 'Argument too small' )
        ifail = 1
c special arguments
      else if( abs( x + dble( m ) ) .lt. xsmall )then
        Gammaf = gbig
        if( ifail .eq. 0 )call crash( 'GAMMAF', 'Argument -ve int' )
        ifail = 1
      else
c factorize as much as possible
        m = x
        if( x .gt. 0.d0 )then
          t = x - dble( m )
          m = m - 1
          g = 1.d0
          if( m .lt. 0 )then
            g = g / x
          else if( m .gt. 0 )then
            do i = 1, m
              g = ( x - dble( i ) ) * g
            end do
          end if
        else
          t = x - dble( m - 1 )
          m = 1 - m
          g = x
          do i = 1, m
            g = ( dble( i ) + x ) * g
          end do
          g = 1.d0 / g
        end if
c Chebyshev expansion evaluated as y(t) for irreducible factor
        t = 2.d0 * t - 1.d0
        y = ( ( ( ( ( ( ( ( ( ( ( ( ( ( ( -1.46381209600000000d-11
     +   * t + 4.26560716800000000d-11 ) * t - 4.01499750400000000d-11 )
     +   * t + 1.27679856640000000d-10 ) * t - 6.13513953280000000d-10 )
     +   * t + 1.82243164160000000d-9 ) * t - 5.11961333760000000d-9 )
     +   * t + 1.53835215257600000d-8 ) * t - 4.64774927155200000d-8 )
     +   * t + 1.39383522590720000d-7 ) * t - 4.17808776355840000d-7 )
     +   * t + 1.25281466396672000d-6 ) * t - 3.75499034136576000d-6 )
     +   * t + 1.12524642975590400d-5 ) * t - 3.36375833240268800d-5 )
     +   * t + 1.00928148823365120d-4 ) * t - 2.96890121633200000d-4
        y = ( ( ( ( ( ( y * t + 9.15785997288933120d-4 ) * t -
     +       2.42259538436268176d-3 )
     +   * t + 9.04033494028101968d-3 ) * t - 1.34118505705967765d-2 )
     +   * t + 1.03703363422075456d-1 ) * t + 1.61691987244425092d-2 )
     +   * t + 8.86226925452758013d-1
c
        Gammaf = y * g
        ifail = 0
      end if
      return
      end
