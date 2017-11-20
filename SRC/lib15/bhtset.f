      subroutine bhtset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to read input from the runXXX.dat file for Barnes-Hut tree code
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
      include 'inc/pi.f'
c
      if( .not. bht )call crash( 'BHTSET', 'Wrong method' )
c read softening length and type
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'SOFT' )then
        print *, 'Next card is'
        print *, line
        call crash( 'BHTSET', 'softening length card not found' )
      end if
      read( line( 11:80 ), *, iostat = i )softl, tsoft
c make cubic spline softening the default
      if( i .ne. 0 )then
        tsoft = 2
        read( line( 11:80 ), *, err = 1 )softl
      end if
      if( softl .le. 0. )call crash( 'BHTSET',
     +                       'Softening length must be greater than 0' )
      if( tsoft .ne. 2 )call crash( 'BHTSET',
     +                            'Cubic spline softening only so far' )
c read opening angle, or tolerance
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'TOLE' )then
        print *, 'Next card is'
        print *, line
        call crash( 'BHTSET', 'tolerance found' )
      end if
      read( line( 11:80 ), *, err = 1 )bhtol
      if( bhtol .le. 0. )call crash( 'BHTSET',
     +                              'Tolerance must be greater than 0' )
      return
c
    1 print *, 'Problem with data card:'
      print *, line( 1:60 )
      call crash( 'BHTSET', 'Error reading data' )
      stop
      end
