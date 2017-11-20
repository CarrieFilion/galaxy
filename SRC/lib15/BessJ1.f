      real*8 function BessJ1( x, ifail )
      implicit none
c returns the value of the Bessel fn J_1 for the given argument x
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
      data xbig, xvsmal / 1.D+16, 3.2D-9 /
C
      t = abs( x )
      if( t .gt. xbig )then
        if( ifail .eq. 0 )call crash( 'BessJ1', 'Argument too large' )
        ifail = 1
        BessJ1 = sqrt( 2.d0 / ( pi * t ) )
      else
        ifail = 0
c small x
        if( t .le. 8.d0 )then
          y = 4.d0
c very small x
          if( t .gt. xvsmal )then
            t = 3.125d - 2 * t * t - 1.d0
            t2 = 2.d0 * t
c expansion evaluated as y(t)  - precision 17e.18
            a = 2.95000000000000000d-18
            b = t2 * a - 1.95540000000000000d-16
            c = t2 * b - a + 1.13857200000000000d-14
            a = t2 * c - b - 5.77740420000000000d-13
            b = t2 * a - c + 2.52812366400000000d-11
            c = t2 * b - a - 9.42421298160000000d-10
            a = t2 * c - b + 2.94970700727800000d-8
            b = t2 * a - c - 7.61758780540030000d-7
            c = t2 * b - a + 1.58870192399321300d-5
            a = t2 * c - b - 2.60444389348580680d-4
            b = t2 * a - c + 3.24027018268385747d-3
            c = t2 * b - a - 2.91755248061542077d-2
            a = t2 * c - b + 1.77709117239728283d-1
            b = t2 * a - c - 6.61443934134543253d-1
            c = t2 * b - a + 1.28799409885767762d+0
            a = t2 * c - b - 1.19180116054121687d+0
            y = t * a - c + 6.48358770605264921d-1
          end if
          BessJ1 = y * x * 0.125d0
        else
c large x
          g = t - .75d0 * pi
          y = sign( sqrt( 2.d0 / ( pi * t ) ), x )
          cx = cos( g ) * y
          sx = -sin( g ) * y * 8.d0 / t
          t = 128.d0 / ( t * t ) - 1.d0
c expansion evaluated as y(t)  - precision 17e.18
          y = ( ( ( ( ( ( ( ( ( ( ( ( ( 1.24928000000000000d-14 )
     +   * t - 3.54508800000000000d-14 ) * t + 6.86387200000000000d-14 )
     +   * t - 2.63372800000000000d-13 ) * t + 1.14622464000000000d-12 )
     +   * t - 5.19113984000000000d-12 ) * t + 2.72040729600000000d-11 )
     +   * t - 1.71310758400000000d-10 ) * t + 1.36030580128000000d-9 )
     +   * t - 1.47085129889600000d-8 ) * t + 2.45367662227560000d-7 )
     +   * t - 7.95969469843846000d-6 ) * t +
     +         8.98804941670557880d-4 ) * t + 1.00090702627808217d+0
c expansion evaluated as g(t)  - precision 17e.18
          g = ( ( ( ( ( ( ( ( ( ( ( ( ( ( ( -2.29376000000000000d-15 )
     +   * t + 5.32480000000000000d-15 ) * t - 4.83328000000000000d-15 )
     +   * t + 1.75718400000000000d-14 ) * t - 7.43936000000000000d-14 )
     +   * t + 2.50224640000000000d-13 ) * t - 9.23228160000000000d-13 )
     +   * t + 3.87554432000000000d-12 ) * t - 1.85248000000000000d-11 )
     +   * t + 1.03950931840000000d-10 ) * t - 7.15110504800000000d-10 )
     +   * t + 6.42013250344000000d-9 ) * t - 8.29196070929200000d-8 )
     +   * t + 1.82120185123076000d-6 ) * t -
     +         9.62145882205441600d-5 ) * t + 4.67768740274489776d-2
          BessJ1 = y * cx + g * sx
        end if
      end if
      return
      end
