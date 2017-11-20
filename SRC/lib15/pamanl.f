      subroutine pamanl( lw, w )
c  Copyright (C) 2015, Jerry Sellwood
c
c performs a Fourier of the mass distribution on the axisymmetric polar grid
c   The grid must be doubled to remove the effects of the periodic images
c   The results overwrite and extend the input array
c
c this version uses single precision FFTPAK routines
      use aarrays
      implicit none
c
c calling arguments
      integer lw
      real w( lw )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local allocatable array
      integer ltrig
      real, allocatable :: trig(:)
c
c local variables
      integer i, j, k, l, mtrig
c
c rearrange and double vertical extent of mass array
      l = 0
      do j = 1, nr( jgrid )
        k = j
        do i = 1, ngz
          l = l + 1
          w( l ) = grdmss( k, 1 )
          w( l + ngz ) = 0
          k = k + nr( jgrid )
        end do
        l = l + ngz
      end do
c initialize
      mtrig = 2 * ngz
      ltrig = 2 * mtrig + 15
      allocate ( trig( ltrig ) )
      call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
c Fourier analysis of each radius in turn
      j = 1
      do i = 1, nr( jgrid )
        call sfftf1( mtrig, w( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c return allocated space
      deallocate ( trig )
      return
      end
