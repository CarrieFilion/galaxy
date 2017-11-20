      real*8 function bk1( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns Bessel K_1( x ) in the range  0 : x : infinity
c   formulae from Abramowitz and Stegun p378
c
c calling argument
      real*8 x
c
c external
      real*8 bi1
c
c local variables
      integer i
      real*8 t, t2, t4, t8
      include 'inc/pi.f'
c
      if( x .lt. 2.0d0 )then
c small x
        t = x / 2.0
        t2 = t * t
        t4 = t2 * t2
        t8 = t4 * t4
        bk1 = 1.  +  0.15443144 * t2
     +            -  0.67278579 * t4
     +            -  0.18156897 * t2 * t4
     +            -  0.01919402 * t8
     +            -  0.00110404 * t2 * t8
     +            -  0.00004686 * t4 * t8
        bk1 = ( bk1 + x * log( t ) * bi1( x ) ) / x
      else if( x .lt. 10.d0 )then
c moderate x
        t = 2.0 / x
        t2 = t * t
        t4 = t2 * t2
        bk1  =  1.25331414
     +       +  0.23498619 * t
     +       -  0.03655620 * t2
     +       +  0.01504268 * t * t2
     +       -  0.00780353 * t4
     +       +  0.00325614 * t * t4
     +       -  0.00068245 * t2 * t4
        bk1 = bk1 / ( sqrt( x ) * exp( x ) )
      else
c asymptotic expansion - eq 9.7.2
        bk1 = 1
        t = 1
        i = 0
        do while ( abs( t ) .gt. 1.d-8 )
          i = i + 1
          t = t * real( 4 - ( 2 * i - 1 )**2 ) / ( real( 8 * i ) * x )
          bk1 = bk1 + t
        end do
        bk1 = bk1 * sqrt( pi / ( 2.d0 * x ) ) * exp( -x )
      end if
      return
      end
