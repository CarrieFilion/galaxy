      subroutine c2manl( l, msstr )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to perform the double Fourier analysis of the mass array for the
c   2-D Cartesian grid.  The Fourier coefficients are stored upwards from
c   msh( ibase ) in an area four times the size of the mesh
c
c this version uses single precision FFTPAK routines
      use aarrays
      implicit none
c
c calling arguments
      integer l
      real msstr( l )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local allocatable arrays
      real, allocatable :: trig(:)
c
      real, allocatable :: w(:)
c
c local variables
      integer i, ix, iy, j, m, mshb, mtrig
c
c allocate space
      mshb = 4 * mesh( jgrid )
      allocate ( w( mshb ) )
c arrange input data and add padding for FFT in y
      j = 0
      do ix = 1, ngx
        m = ix
        do iy = 1, ngy
          j = j + 1
          w( j ) = grdmss( m, 1 )
          w( j + ngy ) = 0
          m = m + ngx
        end do
        j = j + ngy
      end do
c initialize if the vector length is new
      i = 4 * max( ngx, ngy ) + 15
      allocate ( trig( i ) )
      mtrig = 2 * ngy
      call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
c Fourier analysis in y
      j = 1
      do i = 1, ngx
        call sfftf1( mtrig, w( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c re-arrange data and add more padding for FFT in x
      j = 0
      do iy = 1, mtrig
        m = iy
        do ix = 1, ngx
          j = j + 1
          msstr( j ) = w( m )
          msstr( j + ngx ) = 0
          m = m + 2 * ngx
        end do
        j = j + ngx
      end do
c initialize if the vector length is new
      if( 2 * ngx .ne. mtrig )then
        mtrig = 2 * ngx
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
      end if
c Fourier analysis in x
      j = 1
      do i = 1, 2 * ngy
        call sfftf1( mtrig, msstr( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
      deallocate ( w )
      deallocate ( trig )
      return
      end
