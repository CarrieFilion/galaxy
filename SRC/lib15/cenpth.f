      subroutine cenpth
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c fits 2nd order polynomials to the paths of the grid centres from the past n
c   centroid values.  The fits are used to create an estimate of the expected
c   motion of the centroid over the next few substeps.  The curve fitted is
c
c       x( t ) = a t^2 + b t + c
c
c   in each coordinate, for each grid centroid.
c
c the chi^2 minimization gives decreasing weight to points more distant
c   in time from the present
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c local arrays - npast is the number past positions used
      real*8 t( npast ), w( npast )
c matrix, vectors and work areas for fit
      real*8 a( 3, 3 ), aa( 3, 3 ), b( 3 )
c
c local variables
      integer i, iuse, j, k, jg
      real*8 x
      save aa, iuse, t, w
      data iuse / 0  /
c
      if( iuse .eq. 0 )then
c set times (in units of recentering interval) and weights
        t( 1 ) = 0
        w( 1 ) = 2
        do i = 2, npast
          t( i ) = real( 1 - i )
          w( i ) = 1. / t( i )**2
        end do
c generate matrix
        do j = 1, 3
          do i = 1, 3
            aa( i, j ) = 0
          end do
        end do
        do i = 1, npast
          aa( 1, 1 ) = aa( 1, 1 ) + w( i )**2 * t( i )**4
          aa( 1, 2 ) = aa( 1, 2 ) + w( i )**2 * t( i )**3
          aa( 1, 3 ) = aa( 1, 3 ) + w( i )**2 * t( i )**2
          aa( 2, 1 ) = aa( 2, 1 ) + w( i )**2 * t( i )**3
          aa( 2, 2 ) = aa( 2, 2 ) + w( i )**2 * t( i )**2
          aa( 2, 3 ) = aa( 2, 3 ) + w( i )**2 * t( i )
          aa( 3, 1 ) = aa( 3, 1 ) + w( i )**2 * t( i )**2
          aa( 3, 2 ) = aa( 3, 2 ) + w( i )**2 * t( i )
          aa( 3, 3 ) = aa( 3, 3 ) + w( i )**2
        end do
c set initial cenfit value - good even on a restart
        iuse = 1
        do jg = 1, ngrid
          do j = 1, 3
            cenfit( 1, j, jg ) = pastp( j, 1, jg )
         end do
        end do
      end if
c work over active grids
      do jg = 1, ngrid
c work over 3 spatial coordinates
        do j = 1, 3
          do i = 1, 3
            b( i ) = 0
          end do
          do i = 1, npast
            x = pastp( j, i, jg )
            b( 1 ) = b( 1 ) + x * w( i )**2 * t( i )**2
            b( 2 ) = b( 2 ) + x * w( i )**2 * t( i )
            b( 3 ) = b( 3 ) + x * w( i )**2
          end do
c copy matrix - this solves the same matrix over and over, at minimal cost
          do k = 1, 3
            do i = 1, 3
              a( i, k ) = aa( i, k )
            end do
          end do
c solve for coefficients and save
          call mat3x3( a, b )
c preserve old predicted centre and predict the next at t=1
          cenfit( 2, j, jg ) = cenfit( 1, j, jg )
          cenfit( 1, j, jg ) = b( 1 ) + b( 2 ) + b( 3 )
        end do
      end do
      return
      end
