      subroutine dr3set
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to read input from the runXXX.dat file for direct-N code
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
      if( .not. dr3d )call crash( 'DR3SET', 'Wrong method' )
c read softening length
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'SOFT' )then
        print *, 'Next card is'
        print *, line
        call crash( 'DR3SET', 'softening length card not found' )
      end if
      read( line( 11:80 ), *, iostat = i )softl, tsoft
c make cubic spline softening the default
      if( i .ne. 0 )then
        tsoft = 2
        read( line( 11:80 ), *, err = 1 )softl
      end if
      if( softl .le. 0. )call crash( 'DR3SET',
     +                       'Softening length must be greater than 0' )
      if( ( tsoft .lt. 1 ) .or. ( tsoft .gt. 3 ) )call crash( 'DR3SET',
     +                              'Softening rule must be 1, 2 or 3' )
c default offset for ptcls file
      kdrct = 0
      heavies = ngrid .gt. 1
      stphev = .not. heavies
      return
c
    1 print *, 'Problem with data card:'
      print *, line( 1:60 )
      call crash( 'DR3SET', 'Error reading data' )
      stop
      end
