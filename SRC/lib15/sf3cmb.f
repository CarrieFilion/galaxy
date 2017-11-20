      subroutine sf3cmb
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to sum masses in different zones on the grid prior to force
c   determination.
c
c In addition it combines the SFP expansion coeffs from different planes and
c   shells for S3D.  It is much more efficient to do this in one sweep
c   than for each particle as it is being processed.
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c local array
      real, allocatable :: coef(:)
c
c local variables
      integer ig, j, k, l, m, n
      real ak, e
c
c loop over methods
      do ig = 1, ngrid
        call switch( ig )
c 3-D SFP method - combine terms from the different planes of the grid
        if( sf3d )then
          if( basset .ne. 'bess' )call crash( 'SF3CMB',
     +                                              'Impossible basis' )
          allocate ( coef( ngz ) )
c work over functions
          do m = -maxm, maxm
            do n = 0, maxn
              ak = n * deltak
c save coefficients in a scratch array
              k = n * ngx + m + maxm + 1
              do j = 1, ngz
                coef( j ) = sfpfld( k, 1 )
                k = k + ngxy
              end do
c work over field planes
              k = n * ngx + m + maxm + 1
              do j = 1, ngz
                sfpfld( k, 1 ) = 0
                sfpfld( k, 2 ) = 0
c sum over source planes and include exponential dilution factors
                do l = 1, ngz
                  e = exp( -ak * abs( j - l ) * dzg )
                  sfpfld( k, 1 ) = sfpfld( k, 1 ) + coef( l ) * e
                  if( l .lt. j )then
                    sfpfld( k, 2 ) = sfpfld( k, 2 ) - coef( l ) * ak * e
                  else if( l .gt. j )then
                    sfpfld( k, 2 ) = sfpfld( k, 2 ) + coef( l ) * ak * e
                  end if
                end do
                k = k + ngxy
              end do
            end do
          end do
          deallocate ( coef )
        end if
      end do
      call switch( 0 )
      return
      end
