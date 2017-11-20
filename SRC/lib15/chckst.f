      subroutine chckst
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c generic routine to set test masses on a grid
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      if( p2d )then
        call p2ckst
      else if( p3d )then
        call p3ckst
      else if( c2d )then
        call c2ckst
      else if( p3a )then
        call packst( .true. )
      else
        call crash( 'CHCKST', 'Unrecognized grid' )
      end if
c combine masses in different zones and possibly from different processors
      mzone = nzones
      call mascmb
      return
      end
