      subroutine sphbjz( lmax )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c sets the first n zeros of spherical Bessel functions j_0 -> j_lmax
c
c As noted on p 444 of Abromowitz & Stegun, the zeros of j_l(x) are the
c   same as the zeros of J_{l+1/2}(x) which are tabulated on p 467
c
c the number of zeros is set by the parameter mjn in the common block
c
c calling argument
      integer lmax
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/sphrbz.f'
c
c external
      real*8 sphbsj
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind, ir, l, l1, n
      real*8 x1, x2, fx, tol
      include 'inc/pi.f'
c
c approximate values - from Abramowitz & Stegun p467
      data ( zeros( n, 1 ), n = 1, 6 ) /
     + 3.141593,  6.283185,  9.424778, 12.566370, 15.707963, 18.849556 /
      data ( zeros( n, 2 ), n = 1, 6 ) /
     + 4.493409,  7.725252, 10.904122, 14.066194, 17.220755, 20.371303 /
      data ( zeros( n, 3 ), n = 1, 6 ) /
     + 5.763459,  9.095011, 12.322941, 15.514603, 18.689036, 21.853874 /
      data ( zeros( n, 4 ), n = 1, 6 ) /
     + 6.987932, 10.417119, 13.698023, 16.923621, 20.121806, 23.304247 /
      data ( zeros( n, 5 ), n = 1, 6 ) /
     + 8.182561, 11.704907, 15.039655, 18.301256, 21.525418, 24.727566 /
      data ( zeros( n, 6 ), n = 1, 6 ) /
     + 9.355812, 12.966530, 16.354710, 19.653152, 22.904551, 26.2 /
      data ( zeros( n, 7 ), n = 1, 6 ) /
     + 10.512835, 14.207392, 17.647975, 20.983463, 24.262768, 27.6 /
      data ( zeros( n, 8 ), n = 1, 6 ) /
     + 11.657032, 15.431289, 18.922999, 22.295348, 25.5, 29. /
c
      if( lmax .gt. mjl - 1 )then
        print *, 'Calling arg', lmax, ' too large, max is', mjl - 1
        call crash( 'SPHBJZ', 'Calling argument too large' )
      end if
c
      do l = 0, lmax
        l1 = l + 1
        do n = 1, mjn
          if( n .gt. 6 )zeros( n, l1 ) = zeros( n - 1, l1 ) + pi
c refine zero
          x1 = zeros( n, l1 ) - .5
          x2 = zeros( n, l1 ) + .5
          tol = 1.d-15
          ir = 0
          ind = 1
          ifail = 1
          do while ( ind .ne. 0 )
            call fndzro( x1, x2, fx, tol, ir, w, ind, ifail )
            fx = sphbsj( x1, l )
          end do
          if( ( ifail .gt. 0 ) .and. ( ifail .lt. 4 ) )then
            print *, 'ifail =', ifail, ' for l, n =', l, n
            call crash( 'SPHBJZ', 'Failed to find zero' )
          end if
          zeros( n, l1 ) = x1
        end do
      end do
      return
      end
