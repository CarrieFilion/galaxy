      real*8 function frdsfc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine for Cauchy principal value integral for radial force from a disc.
c   It avoids trouble very close to the singularity in frdsfn by building a
c   small table of values a short distance on either side and using
c   interpolated values between the centremost pair only
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      common / feildp / rf, zf
      real*8 rf, zf
c
c external - raw radial force function
      real*8 aitken2, frdsfn
c
c local arrays
      integer n
      parameter ( n = 10 )
      real*8 x( n ), y( n )
c
c local variables
      integer i
      real*8 oldrf, r1, r2, zero
      parameter ( zero = 0 )
      save oldrf, r1, r2
c
      if( rf .ne. oldrf )then
        oldrf = max( rf, 4.5d-3 * rmax )
        do i = 1, n
          x( i ) = oldrf + ( dble( i ) - 5.5 ) * 1.d-3 * rmax
          x( i ) = max( x( i ), zero )
          y( i ) = ( x( i ) - rf ) * frdsfn( x( i ) )
        end do
        oldrf = rf
        r1 = rf - 0.5d-3 * rmax
        r2 = rf + 0.5d-3 * rmax
      end if
c
      if( ( r .gt. r1 ) .and. ( r .lt. r2 ) )then
        frdsfc = aitken2( x, y, n, r )
      else
        frdsfc = ( r - rf ) * frdsfn( r )
      end if
      return
      end
