      real*8 function lgsigmi( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the log of the disk surface density, computed from integrating
c    the DF over all allowed velocities, at the radius r
c
c calling argument
      real*8 r
c
c external
      real*8 gsigmi
c
c local variables
      real*8 r2, sigma
c
      r2 = r
      sigma = gsigmi( r2 )
      lgsigmi = -10.
      if( sigma .gt. 0. )lgsigmi = log( sigma )
      return
      end
