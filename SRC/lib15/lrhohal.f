      real*8 function lrhohal( lr )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the logarithm (base 10) of the halo density at the radius 10**lr
c
c calling argument
      real*8 lr
c
c external
      real*8 rhohal
c
c local variables
      real*8 r, rho
c
      r = 10.**lr
      rho = rhohal( r )
      lrhohal = -10
      if( rho .gt. 0.d0 )lrhohal = log10( rho )
      return
      end
