      subroutine p2fsyn( itype, jr, kr )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to re-synthesize the Fourier transformed force fields or
c   potentials for the 2-D polar grid.  The synthesized data
c   overwrite the input coefficients
c the calling arguments are the type of the field array and the inner
c   and outer radii of the partial grid that is active (could also be
c   the boundaries)
c
c This version uses single precision FFTPAK routines
      use aarrays
      implicit none
c
c calling arguments
      integer itype, jr, kr
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
      integer i, j, k
      real strig
c
c allocate space
      ltrig = 2 * na + 15
      allocate ( trig( ltrig ) )
c initialize trig array
      call sffti1( na, trig( na + 1 ), trig( 2 * na + 1 ) )
c set pointer
      j = ( jr - 1 ) * na + 1
c Fourier synthesis of each grid ring in turn
      do i = jr, kr
        call sfftb1( na, grdfld( j, itype ),
     +                        trig, trig( na + 1 ), trig( 2 * na + 1 ) )
        j = j + na
      end do
c normalize
      strig = 1. / real( na )
      j = ( jr - 1 ) * na + 1
      k = kr * na
      do i = j, k
        grdfld( i, itype ) = strig * grdfld( i, itype )
      end do
c return local workspace
      deallocate ( trig )
      return
      end
