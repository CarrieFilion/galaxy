      real*8 function rhohli( r )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c computes halo density from a double integral over distribution fn
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/moment.f'
c
c externals
      real*8 velint
c
c set power factors for velocity weighting
      iu = 0
      iv = 0
      rhohli = velint( r )
      return
      end
