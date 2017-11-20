      subroutine setspc
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, j, k
c
c space for grids
      do k = 1, ngrid
        call switch( k )
        if( p2d .or. c2d .or. p3a )then
          allocate ( grdmss( mesh( jgrid ), nzones ) )
          allocate ( grdfld( mesh( jgrid ), 3 ) )
          if( nzones .gt. 1 )
     +    allocate ( grdmoz( mesh( jgrid ), 2:nzones ) )
        else if( p3d )then
          allocate ( grdmss( mesh( jgrid ), 0:nzones ) )
          allocate ( grdfld( mesh( jgrid ), 4 ) )
          if( nzones .gt. 1 )
     +    allocate ( grdmoz( mesh( jgrid ), 2:nzones ) )
        else if( c3d )then
          allocate ( grdmss( mesh( jgrid ), 1 ) )
          allocate ( grdfld( mesh( jgrid ), 4:4 ) )
          allocate ( c3dfld( 2 * ngxy, 3 ) )
        else if( s3d )then
          i = 1
          if( hybrid )i = 2
          allocate ( s3dmss( mesh( jgrid ), nzones, i ) )
          allocate ( s3dfld( mesh( jgrid ), i ) )
          if( nzones .gt. 1 )
     +    allocate ( s3dmoz( mesh( jgrid ), 2:nzones, i ) )
c space for field methods
        else if( lsfp )then
          allocate ( sfpmss( mesh( jgrid ), nzones ) )
          if( nzones .gt. 1 )
     +    allocate ( sfpmoz( mesh( jgrid ), 2:nzones ) )
          if( sf3d )then
            allocate ( sfpfld( mesh( jgrid ), 2 ) )
          else
            allocate ( sfpfld( mesh( jgrid ), 1 ) )
          end if
c space for direct methods
        else if( bht )then
          bhmt = 5 * nbod
          bhnlev = 30
          allocate ( bhbox( bhnlev ) )
          allocate ( bhcom( 4, bhmt ) )
          allocate ( itup( bhmt ) )
          allocate ( itdown( bhmt / 8 ) )
        else if( ldrct )then
c find total of all particles to be computed by dr3d
          ndrct = 0
          do i = 1, ncmp
            j = igrd( i )
            if( igrid( j ) .eq. 10 )ndrct = ndrct + nsp( i )
          end do
          allocate ( drpt( 15, ndrct ) )
        end if
      end do
c space for particles
      allocate ( ptcls( lptcls ) )
c space for most bound particles, if needed
      if( lbind )allocate ( bindp( 5, nebind, ngrid ) )
      return
      end
