      real*8 function bi1( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns I_1( x ) in the range -3.75 : x : infinity
c   polynomial approximations from Abramowitz and Stegun p378
c
c calling argument
      real*8 x
c
c local variables
      real*8 t, t2, t4, t6, t8
c
      if( x .lt. 3.75 )then
c small x
        t = x / 3.75
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        bi1 = .5  +
     +        0.87890594 * t2  +
     +        0.51498869 * t4  +
     +        0.15084934 * t2 * t4 +
     +        0.02658733 * t8  +
     +        0.00301532 * t2 * t8  +
     +        0.00032411 * t4 * t8
        bi1 = x * bi1
      else
c large x
        t = 3.75 / x
        t2 = t * t
        t4 = t2 * t2
        t6 = t2 * t4
        bi1 = 0.39894228 -
     +        0.03988024 * t -
     +        0.00362018 * t2 +
     +        0.00163801 * t * t2 -
     +        0.01031555 * t4  +
     +        0.02282967 * t * t4 -
     +        0.02895312 * t6 +
     +        0.01787654 * t * t6 -
     +        0.00420059 * t4 * t4
        bi1 = bi1 * exp( x ) / sqrt( x )
      end if
      return
      end
