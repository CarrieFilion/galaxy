      subroutine ftrans
c  Copyright (C) 2016, Jerry Sellwood
      use aarrays
      implicit none
c routine to compute the power spectrum of a time series of data
c   from a simulation
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
c local arrays
      complex, allocatable :: cd( : )
      real, allocatable :: wk( : )
c
c local variables
      integer ip, it, lt, mp, n
      real an
c
      lt = kt - jt + 1
      mp = kp - jp + 1
      allocate ( cd( lt ) )
      allocate ( wk( 4 * lt ) )
c set normalisation constant
      an = lt
      an = 1. / sqrt( an )
      call szffti( lt, wk )
c wk over p
      do ip = 1, mp
c assemble data for fft
        n = ip - mp
        do it = 1, lt
          n = n + mp
          cd( it ) = cmplx( sdata( 1, n ), sdata( 2, n ) )
        end do
        call szfftf( lt, cd, wk )
c replace time series by power spectrum
        n = ip - mp
        do it = 1, lt
          n = n + mp
          sdata( 1, n ) = an * real( cd( it ) )
          sdata( 2, n ) = an * imag( cd( it ) )
        end do
      end do
      return
      end
