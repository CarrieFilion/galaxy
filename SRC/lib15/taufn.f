      real*8 function taufn( r )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c returns reciprocal of radial speed or zero if imaginary
c   factor sqrt( ( rapo - r ) * ( r - rperi ) ) is for Gauss-Jacobi quad rule
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/orbval.f'
c
c external
      real*8 phitot
c
c local variable
      real*8 u2
c
      taufn = 0
      u2 = 2. * ( crrntE - phitot( r ) )
      if( abs( crrntL ) .gt. 0. )then
c function for Gauss-Jacobi quadrature
        if( r .gt. 0.d0 )u2 = u2 - ( crrntL / r )**2
        if( u2 .gt. 0.d0 )taufn =
     +                      sqrt( ( rapo - r ) * ( r - rperi ) / u2 )
      else
c simple function for radial orbit
        if( u2 .gt. 0.d0 )taufn = 1. / sqrt( u2 )
      end if
      return
      end
