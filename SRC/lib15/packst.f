      subroutine packst( first )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to create an empty mass array and place a small number of
c   unit masses at points read in from data cards to be used as a
c   test of the field determination of the axisymmetric polar grid
c   The direct convolution is is performed by PACHK2
c the calling argument allows for the mass array to be re-created
c   from the original positions without having to duplicate the input
c   data
      use aarrays
      implicit none
c
c calling argument
      logical first
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
      integer i, j
c
      if( .not. p3a )call crash( 'PACKST', 'Wrong type of grid' )
c clear mass array
      do i = 1, mesh( jgrid )
        grdmss( i, 1 ) = 0.
      end do
c read in data
      if( first )then
        if( nzones
     +   .gt. 1 )call crash( 'PACKST', 'Multiple zones not programmed' )
c read in data
        call getline( ni, line )
        if( line( 1:9 ) .ne. 'TEST MASS' )call crash( 'PACKST',
     +                                        'Test mass card missing' )
        read( line( 11:40 ), * )nmass
        nmass = min( nmass, 10 )
        if( master )write( no, '( i10, '' test masses'' )' )nmass
      end if
c put in specified masses
      do i = 1, nmass
        if( first )then
          call getline( ni, line )
          read( line( 11:40 ), * )r( i ), z( i )
          r( i ) = max( r( i ), 1 )
          r( i ) = min( r( i ), nr( jgrid ) )
          z( i ) = max( z( i ), 1 )
          z( i ) = min( z( i ), ngz )
          if( master )write( no,
     +                         '( '' Position'', 2i5 )' )r( i ), z( i )
        end if
        j = ( z( i ) - 1 ) * nr( jgrid ) + r( i )
        grdmss( j, 1 ) = 1
      end do
      return
      end
