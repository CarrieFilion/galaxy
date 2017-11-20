      real*8 function rhalfn( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c external function for density of a spherical mass distribution
c
c calling argument
      real*8 r
c
c external
      real*8 vhalo
c
      rhalfn = r * vhalo( r )**2
      return
      end
