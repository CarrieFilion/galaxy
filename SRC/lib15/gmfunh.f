      real*8 function gmfunh( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for numerical estimation of M(r) for a spherical mass distribution
c
c calling argument
      real*8 r
c
c external
      real*8 rhohal
c
      gmfunh = r * r * rhohal( r )
      return
      end
