      real*8 function gdsphf( r )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c integrand of a spherical disk volume density at a given spherical radius
c
c calling argument in model units, not grid units
      real*8 r
c
c external
      real*8 rhdsph
c
c local variable
      include 'inc/pi.f'
c
c the external returns the mean spherical density at the given radius
      gdsphf = 4. * pi * r**2 * rhdsph( r )
      return
      end
