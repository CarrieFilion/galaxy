      real*8 function bik1( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns I_1( x ) * K_1( x ) in the range  0 : x : infinity
c   polynomial approximations from Abramowitz and Stegun p378
c
c calling argument
      real*8 x
c
c externals
      real*8 bi1, bk1
c
c local variables
      real*8 si1, sk1, t, t2, t4, t6, t8
c
      if( ( x .gt. 2. ) .and. ( x .lt. 3.75 ) )then
        bik1 = bi1( x ) * bk1( x )
      else if( x .le. 2. )then
c small x
        t = x / 3.75
c polynomial for i1
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        si1 = 0.5 +
     +        0.87890594 * t2 +
     +        0.51498869 * t4 +
     +        0.15084934 * t2 * t4 +
     +        0.02658733 * t8 +
     +        0.00301532 * t2 * t8 +
     +        0.00032411 * t4 * t8
c polynomial for k1
        t = x / 2.0
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        sk1 = 1. +
     +        0.15443144 * t2 -
     +        0.67278579 * t4 -
     +        0.18156897 * t2 * t4 -
     +        0.01919402 * t8 -
     +        0.00110404 * t2 * t8 -
     +        0.00004686 * t4 * t8
        sk1 = sk1 + x * log( t ) * bi1( x )
        bik1 = si1 * sk1
      else
c large x
        t = 3.75 / x
c polynomial for i1
        t2 = t * t
        t4 = t2 * t2
        t6 = t2 * t4
        si1 = 0.39894228 -
     +        0.03988024 * t -
     +        0.00362018 * t2 +
     +        0.00163801 * t * t2 -
     +        0.01031555 * t4 +
     +        0.02282967 * t * t4 -
     +        0.02895312 * t6 +
     +        0.01787654 * t * t6 -
     +        0.00420059 * t4 * t4
c polynomial for k1
        t = 2.0 / x
        t2 = t * t
        t4 = t2 * t2
        sk1 = 1.25331414 +
     +        0.23498619 * t -
     +        0.03655620 * t2 +
     +        0.01504268 * t * t2 -
     +        0.00780353 * t4 +
     +        0.00325614 * t * t4 -
     +        0.00068245 * t2 * t4
        bik1 = si1 * sk1 / x
      end if
      return
      end
