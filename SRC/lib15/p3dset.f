      subroutine p3dset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to read in grid parameters for 3-D polar grid
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
      if( .not. p3d )call crash( 'P3DSET', 'Wrong grid type' )
c read grid size
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'GRID' )then
        print *, 'Next card is'
        print *, line
        call crash( 'P3DSET', 'grid size card not found' )
      end if
      read( line( 11:80 ), * )nr( jgrid ), na, ngz
c read vertical spacing of grid
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'Z SP' )then
        print *, 'Next card is'
        print *, line
        call crash( 'P3DSET', 'z spacing card not found' )
      end if
      read( line( 11:80 ), * )dzg
c read number of sectoral harmonics
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'HASH' )then
        print *, 'Next card is'
        print *, line
        call crash( 'P3DSET',
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
          call crash( 'P3DSET', 'uoffset card not found' )
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
c interpolation scheme
      if( line( 1:4 ) .eq. 'INTE' )then
        read( line( 11:80 ), * )jmass
        if( ( jmass .lt. 2 ) .or.
     +      ( jmass .gt. 3 ) )call crash( 'P3DSET',
     +                           'Interpolation scheme must be 2 or 3' )
        call getline( ni, line )
      else
        jmass = 2
      end if
c softening length
      if( line( 1:4 ) .ne. 'SOFT' )then
        print *, 'Next card is'
        print *, line
        call crash( 'P3DSET', 'softl card not found' )
      end if
      read( line( 11:80 ), *, iostat = i )softl, tsoft
c make cubic spline softening the default
      if( i .ne. 0 )then
        tsoft = 2
        read( line( 11:80 ), * )softl
      end if
      if( softl .le. 0. )call crash( 'P3DSET',
     +                       'Softening length must be greater than 0' )
c softening length should not be less than vertical grid spacing
      if( softl .lt. dzg )call crash( 'P3DSET',
     + 'Softening length cannot be less than z-spacing of grid planes' )
      return
      end
