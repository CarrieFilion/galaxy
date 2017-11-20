      subroutine c2dset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up grid parameters for 2-D Cartesian grid
c
c Reads data cards from the .dat file which has been opened previously as
c   logical unit ni
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
      integer i, j, k
c
      if( .not. c2d )call crash( 'C2DSET', 'Wrong method' )
c read grid size
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'GRID' )then
        print *, 'Next card is'
        print *, line
        call crash( 'C2DSET', 'grid size card not found' )
      end if
      read( line( 11:80 ), * )ngx, ngy
      if( ngx .ne. ngy )
     +      call crash( 'C2GRNM', 'Unequal grid dimensions has a bug!' )
c ensure that grid dimensions are suitable for FFTs
      i = 0
      do while ( mod( ngx, 2 ) .eq. 0 )
        i = i + 1
        ngx = ngx / 2
      end do
      j = 0
      do while ( mod( ngx, 3 ) .eq. 0 )
        j = j + 1
        ngx = ngx / 3
      end do
      k = 0
      do while ( mod( ngx, 5 ) .eq. 0 )
        k = k + 1
        ngx = ngx / 5
      end do
      if( ngx .ne. 1 )then
        if( master )write( no, '( ''ngx ='', i4, '' x 2**'',' //
     +        ' i1, '' x 3**'', i1, '' x 5**'', i1 )' )ngx, i, j, k
        call crash( 'C2DSET', 'ngx must be a multiple of 2, 3 & 5' )
      end if
      ngx = 2**i * 3**j * 5**k
      i = 0
      do while ( mod( ngy, 2 ) .eq. 0 )
        i = i + 1
        ngy = ngy / 2
      end do
      j = 0
      do while ( mod( ngy, 3 ) .eq. 0 )
        j = j + 1
        ngy = ngy / 3
      end do
      k = 0
      do while ( mod( ngy, 5 ) .eq. 0 )
        k = k + 1
        ngy = ngy / 5
      end do
      if( ngy .ne. 1 )then
        if( master )write( no, '( ''ngy ='', i4, '' x 2**'',' //
     +        ' i1, '' x 3**'', i1, '' x 5**'', i1 )' )ngy, i, j, k
        call crash( 'C2DSET', 'ngy must be a multiple of 2, 3 & 5' )
      end if
      ngy = 2**i * 3**j * 5**k
c softening length
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'SOFT' )then
        print *, 'Next card is'
        print *, line
        call crash( 'C2DSET', 'softl card not found' )
      end if
      read( line( 11:80 ), *, iostat = i )softl, tsoft
      if( i .ne. 0 )then
c default softening type is Plummer in 2D
        tsoft = 1
        read( line( 11:80 ), * )softl
      end if
      if( softl .le. 0. )
     + call crash( 'C2DSET', 'Softening length must be greater than 0' )
      if( tsoft .ne. 1
     +            )call crash( 'C2DSET', 'Plummer softening rule only' )
      return
      end
