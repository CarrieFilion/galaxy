      real*8 function dfkaln( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Distribution function of the form e**(m-1)*g(x)
c   as Kalnajs (1976, ApJ v205, p751)
c
c calling arguments
      real*8 E, Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 actj1, Ej1j2, gofx
c
c local variables
      real*8 aj1, dfm, dfret, E0, x
c
      dfm = dfcns( 1, icmp )
      if( retract )then
c retrograde star rule - as in IAU 77, p 119 and e-mail 20 Mar 1992
        aj1 = actj1( E, Lz )
        E0 = Ej1j2( aj1 + abs( Lz ), 0.d0 )
        dfret = .5d0 * ( E0 / Phi0 )**( dfm - 1. ) * gofx( 0.d0 )
        if( Lz .gt. 0.d0 )then
c evaluate x
          x = -Lz * sqrt( -2. * E ) / ( Phi0 * rstar )
          dfkaln = ( E / Phi0 )**( dfm - 1. ) * gofx( x )
          dfkaln = dfkaln - dfret
        else
          dfkaln = dfret
        end if
      else
c all stars direct - for the purposes of this routine
        x = -Lz * sqrt( -2. * E ) / ( Phi0 * rstar )
        dfkaln = ( E / Phi0 )**( dfm - 1. ) * gofx( x )
      end if
      return
      end
