      subroutine c2ckst
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c sets up an empty mass array with a few sparsely placed unit masses
c   in readiness for a cross check of the fields calculated by FFTs
c   on the 2-D Cartesian grid
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
      include 'inc/source.f'
c
c local variables
      character*80 line
      integer i, k
c
      if( .not. c2d )call crash( 'C2CKST', 'Wrong type of grid' )
c clear mass array
      do i = 1, mesh( jgrid )
        grdmss( i, 1 ) = 0.
      end do
      if( nzones
     +   .gt. 1 )call crash( 'C2CKST', 'Multiple zones not programmed' )
c read in data
      call getline( ni, line )
      if( line( 1:9 ) .ne. 'TEST MASS' )call crash( 'C2CKST',
     +                                        'Test mass card missing' )
      read( line( 11:40 ), * )nmass
      nmass = min( nmass, 10 )
      if( master )write( no, '( i10, '' test masses'' )' )nmass
c put in specified masses
      do i = 1, nmass
        call getline( ni, line )
        read( line( 11:40 ), * )r( i ), th( i )
        r( i ) = max( r( i ), 1 )
        r( i ) = min( r( i ), ngx )
        th( i ) = max( th( i ), 1 )
        th( i ) = min( th( i ), ngy )
        if( master )write( no, '( '' Position'', 2i5 )' )r( i ), th( i )
        k = ( th( i ) - 1 ) * ngx + r( i )
        grdmss( k, 1 ) = grdmss( k, 1 ) + 1
      end do
      return
      end
