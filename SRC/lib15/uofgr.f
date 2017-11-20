      real function uofgr( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Determines the grid coordinate u for the input cylindrical radius in
c   grid units.
c Unlike the very similar function UOFR, the functional relation between
c   the returned value and input argument can depend on the input value.
c   ie this routine allows for a different spacing of grid points inside
c   what used to be the "hole"
c
c calling argument
      real r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      if( p2d .or. p3d .or. sf2d )then
        if( stdpol )then
          uofgr = log( 1. + r ) / alpha
        else
c various old rules
          if( logrid )then
            if( ( r .gt. 1. ) .or. hole )then
              uofgr = log( r ) / alpha + uoffs
            else
              uofgr = r * uoffs
            end if
          else
            uofgr = r**ibeta - uoffs
          end if
        end if
      else if( p3a )then
        uofgr = r**( 1. / beta )
      else if( s3d )then
        uofgr = log10( r + 1 ) / s3dex
      else if( c2d .or. c3d )then
        uofgr = r
      else
        call crash( 'UOFGR', 'Unrecognised code' )
      end if
      return
      end
