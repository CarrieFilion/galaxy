      real*8 function akcrit( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the critical wavenumber = 2pi/lambda_crit (Toomre 1964) for
c   gravitational instability in a rotationally supported disk
c
c calling argument
      real*8 r
c
c externals
      real*8 akappa, gsigma
c
c local variables
      real*8 ak, s
      include 'inc/pi.f'
c
      s = gsigma( r )
      if( s .gt. 0. )then
        ak = akappa( r )
        akcrit = ak * ak / ( 2. * pi * s )
      else
        akcrit = 1.d10
      end if
      return
      end
