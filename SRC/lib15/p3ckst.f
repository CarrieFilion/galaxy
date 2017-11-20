      subroutine p3ckst
c  Copyright (C) 2015, Jerry Sellwood
c
c sets up an empty mass array with a few sparsely placed unit masses
c   in readiness for a cross check of the fields calculated by FFTs
c   on the 3-D polar grid
      use aarrays
      implicit none
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
      character line*80
      integer i, j, k
c
      if( .not. p3d )call crash( 'P3CKST', 'Wrong type of grid' )
c
      istep = 0
c clear mass arrays
      call switch( 1 )
      do izone = 1, nzones
        do i = 1, mesh( jgrid )
          grdmss( i, izone ) = 0.
        end do
        mstep( izone ) = istep
        jrad( izone ) = nr( jgrid )
        krad( izone ) = 0
        nstep( izone ) = 2**( izone - 1 )
        if( izone .gt. 1 )lstep( 2, izone ) = istep
      end do
c read in data
      call getline( ni, line )
      if( line( 1:9 ) .ne. 'TEST MASS' )call crash( 'P2CKST',
     +                                        'Test mass card missing' )
      read( line( 11:40 ), * )nmass
      nmass = min( nmass, 10 )
      if( master )write( no, '( i10, '' test masses'' )' )nmass
      if( master .and. ( nmass .lt. numprocs ) )print *,
     +          'Fewer test masses than nodes, some grids will be empty'
c put in specified masses
      do i = 1, nmass
        call getline( ni, line )
        read( line( 11:40 ), * )r( i ), th( i ), z( i )
        r( i ) = max( r( i ), 1 )
        r( i ) = min( r( i ), nr( jgrid ) )
        th( i ) = max( th( i ), 1 )
        th( i ) = min( th( i ), na )
        z( i ) = max( z( i ), 1 )
        z( i ) = min( z( i ), ngz )
        izone = r( i ) / 5 + 1
        izone = min( izone, nzones )
        if( master )write( no, '( '' Position'', 3i5, '' zone'', i3 )' )
     +       r( i ), th( i ), z( i ), izone
c assign masses to different nodes
        k = 0
        if( parallel )k = mod( i - 1, numprocs )
        if( k .eq. myid )then
          j = ( ( r( i ) - 1 ) * ngz + z( i ) - 1 ) * na + th( i )
          grdmss( j, izone ) = 1.
c set boundaries of active mass
          jrad( izone ) = min( jrad( izone ), r( i ) )
          krad( izone ) = max( krad( izone ), r( i ) )
        end if
      end do
      return
      end
