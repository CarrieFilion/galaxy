      real*8 function bsjfun( m, x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns J_m(x) - algorithm based on Numerical Recipes p 175
c
c calling arguments
      integer m
      real*8 x
c
c externals
      real*8 BessJ0, BessJ1
c
c parameters
      integer iacc
      real*8 bigno, bigni
      parameter ( iacc = 50, bigno = 1.e10, bigni = 1.e-10 )
c
c local variables
      integer ifail, j, jsum, n
      real*8 ax, bj, bjm, bjp, sum, tox
c
c check for sensible arguments
      if( m .lt. 0 )then
        print *, m, x
        call crash( 'BSJFUN', 'Impossible calling arguments' )
      end if
      ifail = 0
      ax = abs( x )
c special arguments
      if( x .eq. 0.d0 )then
        if( m .eq. 0 )then
          bsjfun = BessJ0( ax, ifail )
        else
          bsjfun = 0
        end if
      else if( m .eq. 0 )then
        bsjfun = BessJ0( ax, ifail )
      else if( m .eq. 1 )then
        bsjfun = BessJ1( ax, ifail )
c use upward recurrence relation where it is stable
      else if( ax .gt. real( m ) )then
        tox = 2. / ax
        bjm = BessJ0( ax, ifail )
        bj = BessJ1( ax, ifail )
        do j = 1, m - 1
          bjp = j * tox * bj - bjm
          bjm = bj
          bj = bjp
        end do
        bsjfun = bj
      else
c downward recurrence from very large m - seems to work!
        tox = 2. / ax
        n = 2 * ( ( m + int( sqrt( real( iacc * m ) ) ) ) / 2 )
        bsjfun = 0
        jsum = 0
        sum = 0
        bjp = 0
        bj = 1
        do j = n, 1, -1
          bjm = j * tox * bj - bjp
          bjp = bj
          bj = bjm
          if( abs( bj ) .gt. bigno )then
            bj = bj * bigni
            bjp = bjp * bigni
            bsjfun = bsjfun * bigni
            sum = sum * bigni
          end if
          if( jsum .ne. 0 )sum = sum + bj
          jsum = 1 - jsum
          if( j .eq. m )bsjfun = bjp
        end do
        sum = 2 * sum - bj
        bsjfun = bsjfun / sum
      end if
      if( ( x .lt. 0. ) .and. ( mod( m, 2 ) .eq. 1 ) )bsjfun = -bsjfun
      return
      end
