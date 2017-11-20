      real*8 function errfn( x, ifail )
      implicit none
c returns the valus of the error function for the given input argument
c   follows the algorithm in the NAG source code
c   the 2nd calling argument, ifail, is redundant
c
c calling arguments - ifail ignored, but included to be consistent with NAG
      integer ifail
      real*8 x
c
c local arrays
      integer nm, ns
      parameter ( nm = 17, ns = 18 )
      real*8 m( nm ), s( ns )
c
c local variables
      integer j
      logical more
      real*8 ax, b, c, d, x2
      include 'inc/pi.f'
c
      data m /  1.4831105640848036d0, -3.010710733865950d-1,
     +          6.89948306898316d-2,  -1.39162712647222d-2,
     +          2.4207995224335d-3,   -3.658639685849d-4,
     +          4.86209844323d-5,     -5.7492565580d-6,
     +          6.113243578d-7,       -5.89910153d-8,
     +          5.2070091d-9,         -4.232976d-10,
     +          3.18811d-11,          -2.2361d-12,
     +          1.467d-13,            -9.0d-15,
     +          5.0d-16 /
      data s /  1.9449071068178803d0, 4.20186582324414d-2,
     +         -1.86866103976769d-2,  5.1281061839107d-3,
     +         -1.0683107461726d-3,   1.744737872522d-4,
     +         -2.15642065714d-5,     1.7282657974d-6,
     +         -2.00479241d-8,       -1.64782105d-8,
     +          2.0008475d-9,         2.57716d-11,
     +         -3.06343d-11,          1.9158d-12,
     +          3.703d-13,           -5.43d-14,
     +         -4.0d-15,              1.2d-15 /
c
      ax = abs( x )
      if( ax .ge. 6.25d0 )then
c large arguments
        errfn = sign( 1.d0, x )
      else if( ax .le. 2.d0 )then
c moderate arguments
        x2 = x * x - 2.d0
        d = 0.d0
        c = m( nm )
        j = nm - 1
        more = .true.
        do while ( more )
          b = x2 * c - d + m( j )
          if( j .gt. 1 )then
            d = c
            c = b
            j = j - 1
          else
            more = .false.
          end if
        end do
        errfn = 0.5d0 * ( b - d ) * x
      else
c small arguments
        x2 = 2.d0 - 20.d0 / ( ax + 3.d0 )
        d = 0.d0
        c = s( ns )
        j = ns - 1
        more = .true.
        do while ( more )
          b = x2 * c - d + s( j )
          if( j .gt. 1 )then
            d = c
            c = b
            j = j - 1
          else
            more = .false.
          end if
        end do
        x2 = 0.5d0 * ( b - d ) * exp( -x * x ) / ( ax * sqrt( pi ) )
        errfn = ( 1.d0 - x2 ) * sign( 1.d0, x )
      end if
      ifail = 0
      return
      end
