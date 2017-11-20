      real*8 function bik0( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns I_0( x ) * K_0( x ) in the range 0 : x : infinity
c   polynomial approximations from Abramowitz and Stegun p378
c
c calling argument
      real*8 x
c
c externals
      real*8 bi0, bk0
c
c local variables
      real*8 si0, sk0, t, t2, t4, t6, t8
c
      if( ( x .gt. 2. ) .and. ( x .lt. 3.75 ) )then
        bik0 = bi0( x ) * bk0( x )
      else if( x .le. 2. )then
c x less than 2
        t = x / 3.75
c polynomial for i0
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        si0 = 1. +
     +        3.5156229 * t2 +
     +        3.0899424 * t4 +
     +        1.2067492 * t2 * t4 +
     +        0.2659732 * t8 +
     +        0.0360768 * t2 * t8 +
     +        0.0045813 * t4 * t8
c polynomial for k0
        t = x / 2.
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        sk0 = -log( t ) * si0
        sk0 = sk0 - 0.57721566 +
     +              0.42278420 * t2 +
     +              0.23069756 * t4 +
     +              0.03488590 * t2 * t4 +
     +              0.00262698 * t8 +
     +              0.00010750 * t2 * t8 +
     +              0.00000740 * t4 * t4
        bik0 = si0 * sk0
      else
c x greater than 3.75
c polynomial for i0
        t = 3.75 / x
        t2 = t * t
        t4 = t2 * t2
        t6 = t2 * t4
        si0 = 0.39894228 +
     +        0.01328592 * t +
     +        0.00225319 * t2 -
     +        0.00157565 * t * t2 +
     +        0.00916281 * t4 -
     +        0.02057706 * t * t4 +
     +        0.02635537 * t6 -
     +        0.01647633 * t * t6 +
     +        0.00392377 * t2 * t6
c polynomial for k0
        t = 2. / x
        t2 = t * t
        t4 = t2 * t2
        sk0 = 1.25331414 -
     +        0.07832358 * t +
     +        0.02189568 * t2 -
     +        0.01062446 * t * t2 +
     +        0.00587872 * t4 -
     +        0.00251540 * t * t4 +
     +        0.00053208 * t2 * t4
        bik0 = si0 * sk0 / x
      end if
      return
      end
