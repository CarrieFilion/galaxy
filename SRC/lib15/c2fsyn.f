      subroutine c2fsyn( itype, l, fldtr )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to perform the double Fourier synthesis of the transformed mesh for
c   the 2-D Cartesian grid.  The results overwrite the input array and are
c   one quarter the size
c
c This version uses single precision FFTPAK routines
c
c calling argument
      integer itype, l
      real fldtr( l )
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
      integer ix, iy, j, m, mtrig
      real scale
c
c allocate space
      j = 2 * mesh( jgrid )
      allocate ( w( j ) )
      j = 4 * max( ngx, ngy ) + 15
      allocate ( trig( j ) )
c initialize
      mtrig = 2 * ngx
      call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
c set pointer
      j = 1
c Fourier synthesis in x first
      do iy = 1, 2 * ngy
        call sfftb1( mtrig, fldtr( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c re-order data for synthesis in y - ignore 2nd half of x-values
      j = 0
      do ix = 1, ngx
        m = ix
        do iy = 1, 2 * ngy
          j = j + 1
          w( j ) = fldtr( m )
          m = m + 2 * ngx
        end do
      end do
c initialize if the vector length is new
      if( mtrig .ne. 2 * ngy )then
        mtrig = 2 * ngy
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
      end if
c set pointer
      j = 1
c Fourier synthesis in y
      do ix = 1, ngx
        call sfftb1( mtrig, w( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c copy results to final area
      m = 0
      scale = 1. / real( 4 * ngx * ngy )
      do iy = 1, ngy
        j = iy
        do ix = 1, ngx
          m = m + 1
          grdfld( m, itype ) = w( j ) * scale
          j = j + 2 * ngy
        end do
      end do
      deallocate ( trig )
      deallocate ( w )
      return
      end
