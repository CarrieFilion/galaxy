      real*8 function hyg2f1( a, b, c, x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c hypergeometric series - defined in terms of Pochammer's symbols
c
c 2f1( a, b; c; x ) = sum_( n = 0 )^( infty )
c                     poch( a )n * poch( b )n / poch( c )n * x**n / n!
c   n.b.  poch( a )n  =  ( a + n - 1 )! / ( a - 1 )!
c
c calling arguments
      real*8 a, b, c, x
c
c local variables
      integer n
      real*8 alarge, ap, bp, cp, term
c
c check for convergence - Matthews & Walker p46
      if( ( abs( x ) .ge. 1.d0 ) .or.
     +    ( abs( x ) .eq. 1.d0 ) .and. ( a + b - c .gt. 0.d0 ) )then
        print 200, a, b, c, x
  200   format( ' args are a =', 1pe12.4, ',  b =', e12.4, ',  c =',
     +        e12.4, ',  x = ', e12.4 )
        call crash( 'HYG2F1', 'non-convergent series in HYG2F1' )
      end if
c use local variables
      ap = a
      bp = b
      cp = c
c initial terms
      term = ap * bp * x / cp
      hyg2f1 = 1.d0 + term
      n = 1
      alarge = abs( term )
c sum power series
      do while ( ( abs( term ) .gt. 1.d-16 * alarge * dble( n ) ) .and.
     +           ( abs( term ) .gt. 1.d-12 * abs( hyg2f1 ) ) )
        n = n + 1
        ap = ap + 1
        bp = bp + 1
        cp = cp + 1
        term = term * ( ap / cp ) * ( bp / dble( n ) ) * x
        alarge = max( alarge, abs( term ) )
        hyg2f1 = hyg2f1 + term
      end do
      return
      end
