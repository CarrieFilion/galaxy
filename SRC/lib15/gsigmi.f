      real*8 function gsigmi( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Computes a double integral over velocities to find surface density
c   - velocity weightings of distfn should not be included
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/moment.f'
c
c external
      real*8 velint
c
      iu = 0
      iv = 0
      gsigmi = velint( r )
      return
      end
