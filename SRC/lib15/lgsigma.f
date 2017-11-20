      real*8 function lgsigma( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the log of the disk surface density at the radius r
c
c calling argument
      real*8 r
c
c external
      real*8 gsigma
c
c local variables
      real*8 sigma
c
      sigma = gsigma( r )
      lgsigma = -10.
      if( sigma .gt. 0. )lgsigma = log( sigma )
      return
      end
