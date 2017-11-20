      real*8 function Omegac( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the circular angular frequency in the current theoretical model
c   at the radius r
c
c calling argument
      real*8 r
c
c externals
      real*8 vcgrad, vcirc
c
c general case
      if( r .eq. 0.d0 )then
        Omegac = vcgrad( 0.d0 )
      else
        Omegac = vcirc( r ) / r
      end if
      return
      end
