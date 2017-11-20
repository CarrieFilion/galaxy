      subroutine p2dset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up grid parameters for 2-D polar grid
c Reads data cards from the .dat file which is opened already as logical
c   unit ni
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
      integer i
c
      if( .not. p2d )call crash( 'P2DSET', 'Wrong type of method' )
c read grid size
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'GRID' )then
        print *, 'Next card is'
        print *, line
        call crash( 'P2DSET', 'grid size card not found' )
      end if
      read( line( 11:80 ), * )nr( jgrid ), na
c read number of sectoral harmonics
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'HASH' )then
        print *, 'Next card is'
        print *, line
        call crash( 'P2DSET',
     +               'highest active sectoral harmonic card not found' )
      end if
      read( line( 11:80 ), * )ng
      ng = ng + 1
c read grid spacing exponent
      call getline( ni, line )
      if( line( 1:4 ) .eq. 'POWE' )then
        read( line( 11:80 ), * )gamma
        if( abs( gamma - 1. ) .gt. 1.e-7 )stdpol = .false.
        call getline( ni, line )
        if( ( .not. stdpol ) .and. ( line( 1:4 ) .ne. 'UOFF' ) )then
          print *, 'Next card is'
          print *, line
          call crash( 'P2DSET', 'uoffset card not found' )
        end if
      else
c set default
        gamma = 1
      end if
c inner edge of grid - default is no inner hole
      if( line( 1:4 ) .eq. 'UOFF' )then
        read( line( 11:80 ), * )i
        uoffs = i
        stdpol = .false.
        call getline( ni, line )
      else
        stdpol = .true.
        uoffs = 0
      end if
c softening length
      if( line( 1:4 ) .ne. 'SOFT' )then
        print *, 'Next card is'
        print *, line
        call crash( 'P2DSET', 'softl card not found' )
      end if
      read( line( 11:80 ), *, iostat = i )softl, tsoft
c make Plummer softening the default
      if( i .ne. 0 )then
        tsoft = 1
        read( line( 11:80 ), * )softl
      end if
      if( softl .le. 0. )call crash( 'P2DSET',
     +                       'Softening length must be greater than 0' )
      return
      end
