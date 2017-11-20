      real function grofu( u )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c determines the cylindrical radius in grid units for the input grid
c   coordinate u
c unlike the very similar function rofu, the functional relation between
c   the returned value and input argument can depend on the input value.
c   ie this routine allows for a different spacing of grid points inside
c   what used to be the "hole"
c
c calling arguments
      real u
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variable
      real x
c
      if( p2d .or. p3d .or. sf2d )then
c obvious rule
        if( stdpol )then
          grofu = exp( alpha * u ) - 1
        else
c original rules
          if( logrid )then
            x = u - uoffs
            if( ( x .ge. 0. ) .or. hole )then
              grofu = exp( alpha * x )
            else
              grofu = u / uoffs
            end if
          else
            grofu = ( u + uoffs )**beta
          end if
        end if
c other grids
      else if( p3a )then
        grofu = u**beta
      else if( s3d )then
        grofu = 10.**( s3dex * u ) - 1
      else if( c2d .or. c3d )then
        grofu = u
      else
        call crash( 'GROFU', 'Unrecognised code' )
      end if
      return
      end
