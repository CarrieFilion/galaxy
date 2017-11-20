      real*8 function lrhohli( lr )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the logarithm (base 10) of the halo density, computed by
c   integration of the DF over all allowed velocities, at the radius 10**lr
c
c calling argument
      real*8 lr
c
c external
      real*8 rhohli
c
c local variables
      real*8 r, rho
c
      r = 10.**lr
      rho = rhohli( r )
      lrhohli = -10
      if( rho .gt. 0.d0 )lrhohli = log10( rho )
      return
      end
