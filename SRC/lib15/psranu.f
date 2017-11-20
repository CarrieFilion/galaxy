      real*8 function psranu( n )
c fills a buffer on the first call with n uniformly distributed, but randomly
c   arranged, fractions between 0 and 1
c   returns one sample, cycling through the buffer on subsequent calls
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling argument
      integer n
c
c external
      real*8 ranuni
c
c local arrays
      integer m
      parameter ( m = 1000 )
      integer ind( m )
      real*8 fr( m )
c
c local variables
      integer i, iuse, j, k, l
      real*8 f
      save fr, iuse, j
c
      data iuse / 0 /
c
c initialize on first call
      if( iuse .ne. n )then
        if( n .gt. m )call crash( 'PSRANU', 'Increase m' )
        do i = 1, n
          ind( i ) = i
        end do
c randomize the array
        k = n
        do i = 1, n
          j = int( dble( k ) * ranuni( 0. ) ) + 1
          l = ind( j )
          ind( j ) = ind( k )
          ind( k ) = l
          k = k - 1
        end do
c set non-random fractions
        f = 1.d0 / real( n )
        do i = 1, n
          j = ind( i )
          fr( j ) = f * ( dble( i ) - .5d0 )
        end do
        j = 0
        iuse = n
      end if
c cycle through buffer
      j = j + 1
      if( j .gt. n )j = 1
      psranu = fr( j )
      return
      end
