      real*8 function gammaP( a, x, ifail )
c from Numerical recipes pp161-163 - ifail is unused, but included for
c   consistency with NAG
      implicit none
c
c calling arguments
      integer ifail
      real*8 a, x
c
c external
      real*8 gammln
c
c local variables
      integer i, itmax
      real*8 an, ap, b, c, d, del, eps, fpmin, gln, h, sum
      parameter ( itmax = 100, eps = 3.d-12, fpmin = 1.e-30 )
c
      if( x .lt. 0. )call crash( 'GAMMAP', 'x is <0' )
      if( a .le. 0. )call crash( 'GAMMAP', 'a is =<0' )
c
      gln = gammln( a )
      if( x .lt. a + 1.d0 )then
c gser from Numerical recipes p162
        ap = a
        sum = 1. / a
        del = sum
        i = 0
        do while ( abs( del ) .gt. abs( sum ) * eps )
          i = i + 1
          if( i .gt. itmax )call crash( 'GAMMAP', 'a too large 1' )
          ap = ap + 1.d0
          del = del * x / ap
          sum = sum + del
        end do
        gammap = sum * exp( a * log( x ) - x - gln )
      else
c gcf from Numerical recipes p162
        b = x + 1.d0 - a
        c = 1.d0 / fpmin
        d = 1.d0 / b
        h = d
        i = 0
        do while ( abs( del - 1.d0 ) .gt. eps )
          i = i + 1
          if( i .gt. itmax )call crash( 'GAMMAP', 'a too large 2' )
          an = -i * ( i - a )
          b = b + 2.
          d = an * d + b
          if( abs( d ) .lt. fpmin )d = fpmin
          c = b + an / c
          if( abs( c ) .lt. fpmin )c = fpmin
          d = 1.d0 / d
          del = d * c
          h = h * del
        end do
        gammap = 1.d0 - exp( a * log( x ) - x - gln ) * h
      end if
      return
      end

      real*8 function gammln( xx )
c from Numerical recipes p157
      implicit none
c
c calling argument
      real*8 xx
c
c local variables
      integer j
      real*8 ser, stp, tmp, x, y, cof( 6 )
c
      data cof, stp / 76.18009172947146d0, -86.50532032941677d0,
     + 24.01409824083091d0, -1.231739572450155d0, .1208650973866179d-2,
     + -.5395239384953d-5, 2.5066282746310005d0 /
c
      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = ( x + 0.5d0 ) * log( tmp ) - tmp
      ser = 1.000000000190015d0
      do j = 1, 6
        y = y + 1.d0
        ser = ser + cof( j ) / y
      end do
      gammln = tmp + log( stp * ser / x )
      return
      end
