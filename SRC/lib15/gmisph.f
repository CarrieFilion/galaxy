      real*8 function gmisph( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns mass as a function of r for spherical isotropic models
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      real*8 aitken2, gmfunh
c
c local arrays
      real*8, allocatable :: absc( : ), wt( : )
c
c local variables
      integer i, j, npts
      real*8 r0, r1
      include 'inc/pi.f'
c
      if( jusedf .eq. 0 )then
        if( master )print *, 'GMISPH: Building a table', nrdd
c rad values may have been set already
        if( ( rmax .gt. 0.d0 ) .and.
     +      ( rtab( nrdd ) .le. rtab( 1 ) ) )then
          do i = 1, nrdd
            rtab( i ) = rmax * real( i - 1 ) / real( nrdd - 1 )
            rtab( i ) = min( rtab( i ), rmax )
          end do
        end if
c initialize mass table
        npts = 64
        allocate ( absc( npts ) )
        allocate ( wt( npts ) )
        mrtab( 1 ) = 0
        do i = 2, nrdd
c non-adaptive Gauss-Legendre quadrature of halo density expression
          r0 = 0
          r1 = rtab( i )
          npts = 32
          call GLqtab( r0, r1, npts, wt, absc )
          mrtab( i ) = 0
          do j = 1, npts
            r0 = absc( j )
            mrtab( i ) = mrtab( i ) + wt( j ) * gmfunh( r0 )
          end do
          mrtab( i ) = 4. * pi * mrtab( i )
        end do
        jusedf = 1
        if( master )print *, 'GMISPH: Table ready'
      end if
c look up value
      if( r .lt. rtab( 2 ) )then
c assume uniform core in innermost bin
        gmisph = r**3 * mrtab( 2 ) / rtab( 2 )**3
      else if( r .gt. rtab( nrdd ) )then
c beyond outer edge
        gmisph = mrtab( nrdd )
      else
c value within tabulated range
        gmisph = aitken2( rtab, mrtab, nrdd, r )
      end if
      return
      end
