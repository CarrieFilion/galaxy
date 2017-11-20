      subroutine p3fsyn( itype, jrf, krf )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to perform the double Fourier synthesis of the field arrays for the
c   3-D polar grid.  The Fourier coefficients overwrite the first half of
c   the input array
c
c this version uses single precision FFTPAK routines
      use aarrays
      implicit none
c
c calling argument
      integer itype, jrf, krf
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
      integer lwork
      real, allocatable :: w(:)
c
c local variables
      integer ia, im, ir, iz, j, m, mtrig, nm, nrf
      real scale
c
c set size of active region
      nrf = krf - jrf + 1
      nm = 0
      do im = 1, ng
        if( .not. lg( im ) )then
          if( im .eq. 1 .or. im .eq. mmax )then
            nm = nm + 1
          else
            nm = nm + 2
          end if
        end if
      end do
c allocate space
      lwork = 2 * nrf * max( 2 * nm, na ) * ngz
      allocate ( w( lwork ) )
      ltrig = 2 * max( na, 2 * ngz ) + 15
      allocate ( trig( ltrig ) )
c initialize if the vector length is new
      if( 2 * ngz .ne. mtrig )then
        mtrig = 2 * ngz
c check space
        j = 2 * mtrig + 15
        if( j .gt. ltrig )call space( ltrig, j, '/ trig /', 'P3FSYN' )
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
      end if
c vertical Fourier synthesis of each vector in turn
      if( itype .lt. 4 )then
        j = 1
        do m = 1, nrf * nm
          call sfftb1( mtrig, grdfld( j, itype ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
          j = j + mtrig
        end do
c re-arrange input data for Fourier synthesis in azimuth
        j = 0
        do ir = 1, nrf
          do iz = 1, ngz
            m = ( ir - 1 ) * nm * mtrig + iz
            do ia = 1, na
              j = j + 1
              w( j ) = 0
              im = ia / 2 + 1
              if( .not. lg( im ) )then
                w( j ) = grdfld( m, itype )
                m = m + mtrig
              end if
            end do
          end do
        end do
      else
c potential coeffs are in grdmss( i, 0 ) array
        j = 1
        do m = 1, nrf * nm
          call sfftb1( mtrig, grdmss( j, 0 ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
          j = j + mtrig
        end do
c re-arrange input data for Fourier synthesis in azimuth
        j = 0
        do ir = 1, nrf
          do iz = 1, ngz
            m = ( ir - 1 ) * nm * mtrig + iz
            do ia = 1, na
              j = j + 1
              w( j ) = 0
              im = ia / 2 + 1
              if( .not. lg( im ) )then
                w( j ) = grdmss( m, 0 )
                m = m + mtrig
              end if
            end do
          end do
        end do
      end if
c initialize if the vector length is new
      if( na .ne. mtrig )then
        mtrig = na
c check space
        j = 2 * mtrig + 15
        if( j .gt. ltrig )call space( ltrig, j, '/ trig /', 'P3FSYN' )
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
      end if
c Fourier synthesis of each ring and plane in turn
      j = 1
      do m = 1, nrf * ngz
        call sfftb1( mtrig, w( j ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        j = j + mtrig
      end do
c copy back results
      scale = 1. / real( 2 * ngz * na )
      m = ngz * na * ( jrf - 1 )
      do j = 1, nrf * ngz * na
        m = m + 1
        grdfld( m, itype ) = scale * w( j )
      end do
      return
      end
