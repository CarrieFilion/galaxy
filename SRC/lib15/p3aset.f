      subroutine p3aset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c reads parameters from the .dat file and set up the axisymmetric polar grid
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c local variables
      character line*80
c
      if( .not. p3a )call crash( 'P3ASET', 'Wrong method' )
c read grid size
      call getline( ni, line )
      read( line, '( 10x, 2i10 )' )nr( jgrid ), ngz
      na = 1
c read grid spacing exponent
      call getline( ni, line )
      read( line, '( 10x, f10.6 )' )gamma
c read shape of grid cells
      call getline( ni, line )
      read( line, '( 10x, f10.6 )' )dzg
      return
      end
