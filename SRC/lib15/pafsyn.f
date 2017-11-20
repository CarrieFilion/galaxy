      subroutine pafsyn( itype, lw, ft )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to re-synthesize the Fourier transformed force fields or
c   potentials for the 2-D polar grid.  The synthesized data
c   overwrite the input coefficients
c
c This version uses single precision FFTPAK routines
      use aarrays
      implicit none
c
c calling arguments
      integer itype, lw
      real ft( lw )
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
      integer ltrig
      real, allocatable :: trig(:)
c
c local variables
      integer i, j, k, l, mtrig
      real strig
c
c initialize
      mtrig = 2 * ngz
      ltrig = 2 * mtrig + 15
      allocate ( trig( ltrig ) )
      call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
      strig = 1. / real( mtrig )
c Fourier synthesis of each radius in turn
      j = 1
      do i = 1, nr( jgrid )
        call sfftb1( mtrig, ft( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c copy results and normalize
      k = 0
      do j = 1, ngz
        l = j
        do i = 1, nr( jgrid )
          k = k + 1
          grdfld( k, itype ) = strig * ft( l )
          l = l + 2 * ngz
        end do
      end do
c return allocated space
      deallocate ( trig )
      return
      end
