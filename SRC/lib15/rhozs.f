      real*8 function rhozs( z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the scaled density as a function of z / z0
c
c calling argument
      real*8 z
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      real*8 ex
      include 'inc/pi.f'
c
      if( ( iztyp( icmp ) .eq. 1 ) .or. ( iztyp( icmp ) .eq. 3 ) )then
c Spitzer sheet - rho(z) = sech^2(z/2z0)/(4z0) & the 4 cancels
        if( abs( z ) .lt. 20.d0 )then
          ex = exp( .5 * z )
          rhozs = 1. / ( ex + 1. / ex )**2
        else
c approximation for large arguments
          rhozs = exp( -abs( z ) )
        end if
      else if( iztyp( icmp ) .eq. 2 )then
c normalized Gaussian distribution
        rhozs = exp( -.5 * z**2 ) / sqrt( 2. * pi )
      else if( iztyp( icmp ) .eq. 4 )then
c rounded exponential distribution
        if( abs( z ) .lt. 4.d0 )then
          rhozs = 1. / ( exp( .5 * abs( z ) ) +
     +                      .2 / exp( 2.5 * abs( z ) ) )**2
        else
c approximation for large arguments
          rhozs = exp( -abs( z ) )
        end if
        rhozs = rhozs / znorm( icmp )
      else if( iztyp( icmp ) .eq. 5 )then
c exponential distribution
        rhozs = .5 * exp( -abs( z ) )
      else
        call crash( 'RHOZS', 'Unknown density profile' )
      end if
      return
      end
