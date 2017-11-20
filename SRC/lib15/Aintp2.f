      subroutine Aintp2( a, b, c, n1, n2, n, x )
      implicit none
c follows the algorithm in the NAG source code for Aitken interpolation
c
c calling arguments
      integer n, n1, n2
      real*8 a( n1 ), b( n1 ), c( n2 ), x
c
c local variables
      integer i, j, k, l, m, nn1
      real*8 dxi, dxk, dyi, dyk
c
      nn1 = n + 1
      do i = 1, nn1
        a( i ) = a( i ) - x
      end do
      m = 2
      k = 1
      l = 0
      do while ( nn1 - k .gt. 0 )
        dxk = a(k)
        dyk = b(k)
        do  i = m, nn1
          j = i + l - 1
          dxi = a( i )
          dyi = ( dyk * dxi - b( i ) * dxk ) / ( dxi - dxk )
          b( i ) = dyi
          c( j ) = dyi
        end do
        l = l + nn1 - m
        k = m
        m = m + 1
      end do
      return
      end
