      subroutine c3dset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to read parameters for 3-D Cartesian grid from .dat file and to
c   set up the input data file for Richard James's Poisson solver
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
c local arrays
      integer ngsize( 3  )
      real*8 wt( 3 )
c
c local variables
      character line*80
      integer i, j, k
      real a
      real*8 den
c
c read in grid size
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'GRID' )then
        print *, 'Record read from unit', ni, ' for grid size'
        print '( a )', line( 1:60 )
        call crash( 'C3DSET', 'Invalid record in .dat file' )
      end if
      read( line( 11:80 ), * )ngx, ngy, ngz
c convert grid dimensions to Richard's input format
      ngsize( 1 ) = ngz
      ngsize( 2 ) = ngy
      ngsize( 3 ) = ngx
      do i = 1, 3
        j = 0
        k = 1
        do while ( k .lt. ngsize( i ) )
          j = j + 1
          k = 2**j + 1
        end do
        if( k .ne. ngsize( i ) )then
          print *, ngsize
          call crash( 'C3DSET', 'Invalid grid dimensions' )
        end if
        ngsize( i ) = j
      end do
c grid cell shapes
      call getline( ni, line )
      if( line( 1:5 ) .ne. 'SHAPE' )then
        print *, 'Record read from unit', ni, ' for grid shape'
        print '( a )', line( 1:60 )
        call crash( 'C3DSET', 'Invalid record in .dat file' )
      end if
      read( line( 11:80 ), * )dh( 3 ), dh( 2 ), dh( 1 )
      dzg = dh( 1 )
c convert shape parameters to Richard's format
      uneq = ( dh( 1 ) .ne. dh( 2 ) ) .or. ( dh( 1 ) .ne. dh( 3 ) )
      a = max( dh( 1 ), dh( 2 ), dh( 3 ) )
      den = 2. * ( dh( 1 )**(-2) + dh( 2 )**(-2) + dh( 3 )**(-2) )
      do i = 1, 3
        wt( i ) = dh( i )**(-2) / den
        dh( i ) = dh( i ) / a
      end do
c mass assignment scheme - CIC is hard wired
      jmass = 2
c scale factor is now 1, and is no longer needed (or used?)
      c3pmf = 1
c set notional softening length - determined from a fit to measured forces
      tsoft = 2
      softl = 1.8
c create Richard's 'steer' file
      call new_unit( i )
      open( i, file = 'steer', form = 'formatted', status = 'unknown' )
      write( i, '( ''MESH      '', 3i10 )' )ngsize
      write( i, '( ''END'' )' )
      close( i )
      return
      end
