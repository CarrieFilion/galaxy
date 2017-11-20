      subroutine checkf
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c generic routine to check solution for field of a few test masses on a grid
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      if( master )then
        if( p2d )then
          call p2chkf
        else if( p3d )then
          call p3chkf
        else if( c2d )then
          call c2chkf
        else if( p3a )then
          call pachkf
        else
          call crash( 'CHECKF', 'Unrecognized grid' )
        end if
      end if
      return
      end
