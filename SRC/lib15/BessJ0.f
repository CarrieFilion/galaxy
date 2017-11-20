      real*8 function BessJ0( x, ifail )
      implicit none
c returns the value of the Bessel fn J_0 for the given argument x
c   follows the algorithm in the NAG source code
c
c calling arguments
      integer ifail
      real*8 x
c
c local variables
      real*8 a, b, c, cx, g, sx, t, t2, xbig, xvsmal, y
      include 'inc/pi.f'
c
      data xbig, xvsmal / 1.d+16, 3.2d-9 /
c
      t = abs( x )
      if( t .gt. xbig )then
        if( ifail .eq. 0 )call crash( 'BessJ0', 'Argument too large' )
        ifail = 1
        BessJ0 = sqrt( 2.d0 / ( pi * t ) )
      else if( t .le. xvsmal )then
        BessJ0 = 1.d0
      else
        ifail = 0
c moderate arguments
        if( t .le. 8.d0 )then
          t = 3.125d-2 * t * t - 1.d0
          t2 = 2.d0 * t
c expansion evaluated as y(t) - precision 17e.18
          a = 1.22200000000000000d-17
          b = t2 * a - 7.58850000000000000d-16
          c = t2 * b - a + 4.12532100000000000d-14
          a = t2 * c - b - 1.94383469000000000d-12
          b = t2 * a - c + 7.84869631400000000d-11
          c = t2 * b - a - 2.67925353056000000d-9
          a = t2 * c - b + 7.60816359241900000d-8
          b = t2 * a - c - 1.76194690776215000d-6
          c = t2 * b - a + 3.24603288210050800d-5
          a = t2 * c - b - 4.60626166206275050d-4
          b = t2 * a - c + 4.81918006946760450d-3
          c = t2 * b - a - 3.48937694114088852d-2
          a = t2 * c - b + 1.58067102332097261d-1
          b = t2 * a - c - 3.70094993872649779d-1
          c = t2 * b - a + 2.65178613203336810d-1
          a = t2 * c - b - 8.72344235285222129d-3
          y = t * a - c + 1.57727971474890120d-1
c
          BessJ0 = y
        else
c large arguments
          g = t - 0.25d0 * pi
          y = sqrt( 2.d0 / ( pi * t ) )
          cx = cos( g ) * y
          sx = -sin( g ) * y * 8.d0 / t
          t = 128.d0 / ( t * t ) - 1.d0
c expansion evaluated as y(t) - precision 17e.18
          y = ( ( ( ( ( ( ( ( ( ( ( ( ( -1.17964800000000000d-14 )
     +   * t + 3.34028800000000000d-14 ) * t - 6.41843200000000000d-14 )
     +   * t + 2.45294080000000000d-13 ) * t - 1.06365696000000000d-12 )
     +   * t + 4.78702080000000000d-12 ) * t - 2.48827276800000000d-11 )
     +   * t + 1.55005642880000000d-10 ) * t - 1.21211819632000000d-9 )
     +   * t + 1.28037614434400000d-8 ) * t - 2.05274481565160000d-7 )
     +   * t + 6.13741608010926000d-6 ) * t
     +       - 5.36367319213004570d-4 ) * t + 9.99457275788251954d-1
c expansion evaluated as g(t) - precision 17e.18
          g = ( ( ( ( ( ( ( ( ( ( ( ( ( ( ( ( -9.83040000000000000d-16 )
     +   * t + 2.12992000000000000d-15 ) * t - 1.14688000000000000d-15 )
     +   * t + 4.75136000000000000d-15 ) * t - 2.27942400000000000d-14 )
     +   * t + 6.95193600000000000d-14 ) * t - 2.28807680000000000d-13 )
     +   * t + 8.59855360000000000d-13 ) * t - 3.59094272000000000d-12 )
     +   * t + 1.70330918400000000d-11 ) * t - 9.47034774400000000d-11 )
     +   * t + 6.43278173280000000d-10 ) * t - 5.66871613024000000d-9 )
     +   * t + 7.10621485930000000d-8 ) * t - 1.47713883264594000d-6 )
     +   * t + 6.83314909934390000d-5 ) * t - 1.55551138795135187d-2
c
          BessJ0 = y * cx + g * sx
        end if
      end if
      return
      end
