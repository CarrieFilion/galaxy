      real*8 function bk0( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns Bessel K_0( x ) in the range    0 : x : infinity
c   formulae from Abramowitz and Stegun p378
c
c calling argument
      real*8 x
c
c external
      real*8 bi0
c
c local variables
      integer i
      real*8 t, t2, t4, t8
      include 'inc/pi.f'
c
      if( x .lt. 2.d0 )then
c small x
        t = x / 2.
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        bk0 = -log( t ) * bi0( x )
        bk0 = bk0  -  0.57721566
     +             +  0.42278420 * t2
     +             +  0.23069756 * t4
     +             +  0.03488590 * t2 * t4
     +             +  0.00262698 * t8
     +             +  0.00010750 * t2 * t8
     +             +  0.00000740 * t4 * t4
      else if( x .lt. 10.d0 )then
c moderate x
        t = 2. / x
        t2 = t * t
        t4 = t2 * t2
        bk0 =    1.25331414
     +        -  0.07832358 * t
     +        +  0.02189568 * t2
     +        -  0.01062446 * t * t2
     +        +  0.00587872 * t4
     +        -  0.00251540 * t * t4
     +        +  0.00053208 * t2 * t4
        bk0 = bk0 / ( sqrt( x ) * exp( x ) )
      else
c asymptotic expansion - eq 9.7.2
        bk0 = 1
        t = 1
        i = 0
        do while ( abs( t ) .gt. 1.d-8 )
          i = i + 1
          t = -t * real( ( 2 * i - 1 )**2 ) / ( real( 8 * i ) * x )
          bk0 = bk0 + t
        end do
        bk0 = bk0 * sqrt( pi / ( 2.d0 * x ) ) * exp( -x )
      end if
      return
      end
