      subroutine smfield( j )
c  Copyright (C) 2017, Jerry Sellwood
      use aarrays
      implicit none
c
c calling argument
      integer j
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      real*8 meanph2, meanrf2, splint2
c
c local variables
      integer i, jcmp, k
      real*8 r, rs
c
      jcmp = icmp
      print *, 'using the grid to compute the field of component', j
c initialize grid parameters etc
      call grdset
      call setrun
      print *, 'done setrun'
      if( .not. lgrd )call crash( 'SMFIELD', 'Not a grid' )
c create or check Greens function file for grid methods
      call greenm( .false. )
      print *, 'done greenm'
c allocate space
      call setspc
c assign mass to grids etc
      phys = .true.
      lprint = .true.
      dist( j ) = .false.
      call smmass
!      if( master )print *, 'done smmass'
c determine gravitational field
      call findf( .true. )
!      if( master )print *, 'done findf'
c create and save tables of potential and central attraction
      k = -1
      call opnfil( k, 'ftb', 'unformatted', 'unknown', 'seq', i )
      if( i .ne. 0 )call crash( 'SMFIELD', 'open ftb file failed' )
      nradad = mradad
      if( .not. allocated( arad ) )allocate ( arad( nradad ) )
      if( .not. allocated( Phin ) )allocate ( Phin( nradad ) )
      if( .not. allocated( frin ) )allocate ( frin( nradad ) )
      do i = 1, nradad
        rs = .999 * rgrid( jgrid ) * real( i - 1 ) / real( nradad - 1 )
        r = rs / lscale
        arad( i ) = r
c these functions include contributions from rigid components
        Phin( i ) = meanph2( rs, 0.d0 ) * gvfac**2
        frin( i ) = meanrf2( rs ) / ( lscale * ts**2 )
      end do
      rewind k
      write( k )nradad, arad, Phin, frin
      close( k )
      icmp = jcmp
c initialize splines
      if( .not. allocated( plamda ) )then
        allocate ( plamda( nradad + 4 ) )
        allocate ( phinc( nradad + 4 ) )
      end if
      if( .not. allocated( flamda ) )then
        allocate ( flamda( nradad + 4 ) )
        allocate ( frinc( nradad + 4 ) )
      end if
      r = arad( nradad / 2 )
      rs = splint2( arad, Phin, nradad, r, plamda, phinc, .true. )
      rs = splint2( arad, frin, nradad, r, flamda, frinc, .true. )
      lgfld = .true.
!      if( master )print *, 'done smfield'
      return
      end
