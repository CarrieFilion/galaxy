      subroutine p2manl( ip, jr, kr )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to perform the Fourier analysis of the mass array for the 2-D
c   polar grid.  The Fourier coefficients overwrite the input array
c the calling paramaters are the base address of the array containing
c   the input data and the inner and outer radii of the partial grid
c   that is active (could also be the boundaries)
c
c This version uses single precision FFTPAK routines
      use aarrays
      implicit none
c
c calling arguments
      integer ip, jr, kr
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
      integer i, j
c
      ltrig = 2 * na + 15
      allocate ( trig( ltrig ) )
c initialize trig array
      call sffti1( na, trig( na + 1 ), trig( 2 * na + 1 ) )
c set base address of active mass region
      j = ( jr - 1 ) * na + 1
c Fourier analysis of each grid ring in turn
      do i = jr, kr
        call sfftf1( na, grdmss( j, ip ),
     +                        trig, trig( na + 1 ), trig( 2 * na + 1 ) )
        j = j + na
      end do
c return local workspace
      deallocate ( trig )
      return
      end
