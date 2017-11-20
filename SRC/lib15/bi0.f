      real*8 function bi0( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns I_0( x ) in the range - 3.75 : x : infinity
c   polynomial approximations from Abramowitz and Stegun p378
c
c calling argument
      real*8 x
c
c local variables
      real*8 t, t2, t4, t6, t8
c
      if( x .le. 3.75 )then
        t = x / 3.75
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        bi0 = 1. +
     +      3.5156229 * t2  +
     +       3.0899424 * t4  +
     +       1.2067492 * t2 * t4 +
     +       0.2659732 * t8  +
     +       0.0360768 * t2 * t8  +
     +       0.0045813 * t4 * t8
      else
c large x
        t = 3.75 / x
        t2 = t * t
        t4 = t2 * t2
        t6 = t2 * t4
        bi0 = 0.39894228 +
     +        0.01328592 * t +
     +        0.00225319 * t2 -
     +        0.00157565 * t * t2 +
     +        0.00916281 * t4 -
     +        0.02057706 * t * t4 +
     +        0.02635537 * t6 -
     +        0.01647633 * t * t6 +
     +        0.00392377 * t2 * t6
         bi0 = bi0 * exp( x ) / sqrt( x )
      end if
      return
      end
