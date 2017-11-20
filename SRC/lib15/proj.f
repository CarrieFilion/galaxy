      real function projz( x, y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c the x,y projection
c
c calling arguments
      real x, y
c
c external
      real proj3d
c
c local variables
      integer idir
      real pos( 3 ), r2
c
      r2 = x**2 + y**2
      if( r2 .ge. 1. )then
        projz = 0
      else
        pos( 2 ) = y
        pos( 3 ) = x
        idir = 1
        projz = proj3d( pos, idir )
      end if
      return
      end

      real function projy( x, z )
c
c calling arguments
      real x, z
c the x,z projection
c
c external
      real proj3d
c
c local variables
      integer idir
      real pos( 3 ), r2
c
      r2 = x**2 + z**2
      if( r2 .ge. 1. )then
        projy = 0
      else
        pos( 1 ) = z
        pos( 3 ) = x
        idir = 2
        projy = proj3d( pos, idir )
      end if
      return
      end

      real function projx( z, y )
c the y,z projection
c
c calling arguments
      real y, z
c
c external
      real proj3d
c
c local variables
      integer idir
      real pos( 3 ), r2
c
      r2 = z**2 + y**2
      if( r2 .ge. 1. )then
        projx = 0
      else
        pos( 1 ) = z
        pos( 2 ) = y
        idir = 3
        projx = proj3d( pos, idir )
      end if
      return
      end
