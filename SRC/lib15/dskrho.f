      real*8 function dskrho( x, y, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Returns the volume density at an arbitrary point.  The position coordinates
c   and returned value are in natural units.
c
c calling arguments in natural units
      real*8 x, y, z
c
c external
      real*8 rhorz
c
c local variables
      real*8 r
c
      r = sqrt( x * x + y * y )
      dskrho = rhorz( r, z )
      return
      end
