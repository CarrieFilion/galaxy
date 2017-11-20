      real*8 function rhorz( r, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c function to return the disk volume density at the point specified in
c   cylindrical polars
c
c calling arguments in model units, not grid units
      real*8 r, z
c
c externals
      real zthick
      real*8 gsigmt, rhozs
c
c local variable
      real*8 z0
c
      z0 = zthick( r )
      rhorz = 0
      if( z0 .gt. 0.d0 )then
        rhorz = gsigmt( r ) * rhozs( z / z0 ) / z0
      end if
      return
      end
