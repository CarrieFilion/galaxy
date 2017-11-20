      real function aitken( x, y, n, xx )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Interpolation by Aitken's method
c        n.b. - works for unevenly spaced abscissae but they must be ordered
c
c calling arguments
      integer n
      real xx, x( n ), y( n )
c
c local arrays
      integer mmax
      parameter ( mmax = 9 )
      integer ix( mmax + 1 )
      real*8 a( mmax + 1 ), b( mmax + 1 ), c( (mmax * (mmax + 1 )) / 2 )
c
c local variables
      integer i, k, m, m1, m2, nmin, nmax
      real*8 x2
c
c find nearest abscissa by bi-section - assumed in ascending order
      nmin = 1
      nmax = n
      do while ( nmin + 1 .lt. nmax )
        i = ( nmin + nmax ) / 2
        if( x( i ) .lt. xx )then
          nmin = i
        else
          nmax = i
        end if
      end do
c choose absolutely nearest abscissa
      if( abs( xx - x( nmin ) ) .lt. abs( xx - x( nmax ) ) )then
        ix( 1 ) = nmin
      else
        ix( 1 ) = nmax
      end if
c choose number of abscissae
      m = min( mmax, n - 1 )
      if( m .eq. 0 )then
        aitken = y( ix( 1 ) )
      else
        m1 = m + 1
c search for next nearest abscissae assuming data are ordered
        nmin = ix( 1 ) - 1
        nmax = ix( 1 ) + 1
        do k = 2, m1
c check for ends of range
          if( ( nmin .gt. 0 ) .and. ( nmax .le. n ) )then
            if( abs( xx - x( nmin ) ) .lt. abs( xx - x( nmax ) ) )then
              ix( k ) = nmin
            else
              ix( k ) = nmax
            end if
          else
            if( nmin .le. 0 )ix( k ) = nmax
            if( nmax .gt. n )ix( k ) = nmin
          end if
          nmax = max( ix( k ) + 1, nmax )
          nmin = min( ix( k ) - 1, nmin )
        end do
c assemble sub-table
        do k = 1, m1
          i = ix( k )
          a( k ) = x( i )
          b( k ) = y( i )
        end do
c interpolation by Aitken's method
        m2 = ( m * m1 ) / 2
        x2 = xx
        call Aintp2( a, b, c, m1, m2, m, x2 )
        aitken = c( m2 )
      end if
      return
      end

      real*8 function aitken2( x, y, n, xx )
c Interpolation by Aitken's method
c        n.b. - works for unevenly spaced abscissae but they must be ordered
c
c calling arguments
      integer n
      real*8 xx, x( n ), y( n )
c
c local arrays
      integer mmax
      parameter ( mmax = 19 )
      integer ix( mmax + 1 )
      real*8 a( mmax + 1 ), b( mmax + 1 ), c( (mmax * (mmax + 1 )) / 2 )
c
c local variables
      integer i, k, m, m1, m2, nmin, nmax
c
c find nearest abscissa by bi-section - assumed in ascending order
      nmin = 1
      nmax = n
      do while ( nmin + 1 .lt. nmax )
        i = ( nmin + nmax ) / 2
        if( x( i ) .lt. xx )then
          nmin = i
        else
          nmax = i
        end if
      end do
c choose absolutely nearest abscissa
      if( abs( xx - x( nmin ) ) .lt. abs( xx - x( nmax ) ) )then
        ix( 1 ) = nmin
      else
        ix( 1 ) = nmax
      end if
c choose number of abscissae
      m = min( mmax, n - 1 )
      if( m .eq. 0 )then
        aitken2 = y( ix( 1 ) )
      else
        m1 = m + 1
c search for next nearest abscissae assuming data are ordered
        nmin = ix( 1 ) - 1
        nmax = ix( 1 ) + 1
        do k = 2, m1
c check for ends of range
          if( ( nmin .gt. 0 ) .and. ( nmax .le. n ) )then
            if( abs( xx - x( nmin ) ) .lt. abs( xx - x( nmax ) ) )then
              ix( k ) = nmin
            else
              ix( k ) = nmax
            end if
          else
            if( nmin .le. 0 )ix( k ) = nmax
            if( nmax .gt. n )ix( k ) = nmin
          end if
          nmax = max( ix( k ) + 1, nmax )
          nmin = min( ix( k ) - 1, nmin )
        end do
c assemble sub-table
        do k = 1, m1
          i = ix( k )
          a( k ) = x( i )
          b( k ) = y( i )
        end do
c interpolation by Aitken's method
        m2 = ( m * m1 ) / 2
        call Aintp2( a, b, c, m1, m2, m, xx )
        aitken2 = c( m2 )
      end if
      return
      end
