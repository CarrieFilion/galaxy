      real*8 function gmfund( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for the numerical estimation of M(R) for a disk mass distribution
c
c calling argument
      real*8 r
c
c external
      real*8 gsigmt
c
c local variable
      include 'inc/pi.f'
c
      gmfund = 2. * pi * r * gsigmt( r )
      return
      end
