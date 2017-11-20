      real*8 function xtoom( r, m )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns Toomre's X parameter at the radius r for the given
c   sectoral harmonic m
c
c calling arguments
      integer m
      real*8 r
c
c externals
      real*8 akappa, gsigma
c
c local variable
      include 'inc/pi.f'
c
      xtoom = r * akappa( r )**2 / ( 2. * pi * real( m ) * gsigma( r ) )
      return
      end
