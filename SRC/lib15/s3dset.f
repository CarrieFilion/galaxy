      subroutine s3dset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to read input from the runXXX.dat file for the PM+SH method
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
      logical firstc
      integer i
      save firstc
      data firstc / .true. /
c
      if( .not. s3d )call crash( 'S3DSET', 'Wrong method' )
c determine grid
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'GRID' )then
        print *, 'Next card is'
        print *, line
        call crash( 'S3DSET', 'grid size card not found' )
      end if
c grid size and outer boundary
      read( line( 11:80 ), * )nr( jgrid )
      if( nr( jgrid ) .gt. mr1
     +               )call space( mr1, nr( jgrid ), 'radius', 's3dset' )
c choose grid outer boundary and order on first call only
      if( firstc )then
        read( line( 11:80 ), *, err = 1 )i, s3rad( nr( jgrid ) )
        call getline( ni, line )
        if( line( 1:4 ) .ne. 'LMAX' )then
          print *, 'Next card is'
          print *, line
          call crash( 'S3DSET', 'lmax card not found' )
        end if
        read( line( 11:80 ), *, err = 1 )s3lmax
        firstc = .false.
      end if
      return
    1 print *, 'problem reading line'
      print '( 1x, a )', line( 1:50 )
      call crash( 'S3DSET', 'Reading input data' )
      stop
      end
