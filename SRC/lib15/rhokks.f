      real*8 function rhokks( r, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c volume density of the Kuz'min-Kutuzov family of oblate spheroids
c
c calling arguments
      real*8 r, z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      real*8 a2, c2, lambda, nu, x
      include 'inc/pi.f'
c
      if( ctype( icmp ) .eq. 'KKSP' )then
        rhokks = 0
        if( sqrt( r**2 + z**2 ) .lt. rmax )then
          call kkcoor( r, z, lambda, nu )
          a2 = rscale( icmp )**2
          c2 = a2 * dfcns( 3, icmp )**2
          x = lambda * nu
c equation (4.3) of Dejonghe & de Zeeuw (1988, ApJ v333 p100)
          rhokks = ( x + a2 * ( lambda + nu + 3. * sqrt( x ) ) ) /
     +                   ( x**1.5 * ( sqrt( lambda ) + sqrt( nu ) )**3 )
          rhokks = c2 * rhokks / ( 4. * pi )
        end if
      else
        call crash( 'RHOKKS', 'Wrong model' )
      end if
      return
      end
