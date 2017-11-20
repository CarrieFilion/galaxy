      subroutine pafndf
c  Copyright (C) 2015, Jerry Sellwood
c
c determines the forces and potentials on the axisymmetric polar grid
c   by FFT in the vertical direction and then by direct convolution
c   in radius.
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
c local allocatable arrays
      real, allocatable :: ft(:)
c
      real, allocatable :: mt(:)
c
c local variables
      integer itype, jty, kty, lw
c
c allocate scratch area for mass transform
      lw = 2 * mesh( jgrid )
      allocate ( mt( lw ) )
c form Fourier transform of mass array
      call pamanl( lw, mt )
c
c allocate scratch area for field transform
      allocate ( ft( lw ) )
c work over types
      if( fixrad )then
        jty = 2
        kty = 2
      else
        jty = 1
        kty = 2
        if( phys )kty = 3
      end if
      do itype = jty, kty
        call blkcpy( mt, ft, lw )
c convolve with Green function
        call paconv( itype, lw, mt, ft )
c re-synthesize
        call pafsyn( itype, lw, ft )
      end do
      deallocate ( ft )
      deallocate ( mt )
      return
      end
