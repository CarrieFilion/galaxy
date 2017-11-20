      real*8 function vmxwel( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Maxwellian velocity distribution as a function of v/sigma
c   Useful for for a 3-D isotropic sphere
c   as eq (2.15) Hernquist (1993, ApJS, v86, p389)
c
c calling argument
      real*8 x
c
c local variable
      include 'inc/pi.f'
c
      vmxwel = 2. * x**2 * exp( -.5 * x**2 ) / sqrt( 2. * pi )
      return
      end
