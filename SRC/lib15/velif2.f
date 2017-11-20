      real*8 function velif2( v )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes velocity weighted distribution function at ( r, u, v )
c
c calling argument
      real*8 v
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
c external
      real*8 distfn
c
c local variables
      real*8 E, Lz
      include 'inc/pi.f'
c
      E = pot + .5 * ( u * u + v * v )
      Lz = rad * v
      velif2 = distfn( E, Lz )
c weight by requested powers of velocity components
      if( disc( icmp ) )then
c disk models
        if( iu .ne. 0 )velif2 = velif2 * u**iu
      else if( sphrod( icmp ) )then
c spheroidal models
        if( iu .eq. 0 )then
c u is not the real radial velocity but the component in the meriodional plane
          velif2 = 2. * pi * u * velif2
        else
c the expectation value of the squared radial velocity is 1/2 u**2 therefore
          velif2 = 2. * pi * u * velif2 * ( 0.5 * u * u )**( iu / 2 )
        end if
      else
c spherical models - the direction of the v vector is isotropic
        velif2 = 2. * pi * v * velif2 * u**iu
      end if
      if( iv .ne. 0 )velif2 = velif2 * v**iv
      return
      end
