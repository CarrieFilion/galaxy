      real*8 function phgenc( rs )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for Cauchy principal value in phgend
c
c It avoids trouble very close to the singular field point by building a
c   small table of values a short distance on either side and using
c   interpolated values between the centremost pair only
c
c calling argument - the radius of the source
      real*8 rs
c
c common block
c
      common / feildp / rf, zf
      real*8 rf, zf
c
c externals - interpolation function and integrand without the ( rs - rf ) term
      real*8 aitken2, phgenf
c
c local arrays
      integer n
      parameter ( n = 10 )
      real*8 x( n ), y( n )
c
c local variables
      integer i
      real*8 oldrf, r1, r2, zero
      parameter ( zero = 0. )
      save oldrf, r1, r2
c
      if( rf .ne. oldrf )then
        do i = 1, n
          x( i ) = rf + ( dble( i ) - 5.5 ) * 1.d-3 * rf
          x( i ) = max( x( i ), zero )
          y( i ) = ( x( i ) - rf ) * phgenf( x( i ) )
        end do
        oldrf = rf
        r1 = rf - 0.5d-3 * rf
        r2 = rf + 0.5d-3 * rf
      end if
c
      if( ( rs .gt. r1 ) .and. ( rs .lt. r2 ) )then
        phgenc = aitken2( x, y, n, rs )
      else
        phgenc = ( rs - rf ) * phgenf( rs )
      end if
      return
      end
