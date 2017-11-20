      real function uofr( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Determines the grid coordinate u for the input cylindrical radius in
c   grid units.
c Unlike the very similar function UOFGR, the functional relation between
c   the returned value and input argument does not depend on the input value.
c   ie this routine makes no allowance for a different spacing of grid points
c   inside what used to be the "hole"
c
c calling arguments
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
      if( p2d .or. sf2d .or. p3d )then
        if( stdpol )then
          uofr = log( 1. + r ) / alpha
        else
c various old rules
          if( logrid )then
            uofr = log( r ) / alpha + uoffs
          else
            uofr = r**ibeta - uoffs
          end if
        end if
      else if( p3a )then
        uofr = r**( 1. / beta )
      else
        call crash( 'UOFR', 'Unrecognised code' )
      end if
      return
      end
