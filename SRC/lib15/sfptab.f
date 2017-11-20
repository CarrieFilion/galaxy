      subroutine sfptab
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c make a look-up table of basis function values for use in SFP methods
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
c externals
      real*8 ajfun, ajder, bsjfun, bsjder, nkaln
c
c local variables
      integer i, j, l, m, n
      real*8 k, r, x
      include 'inc/pi.f'
c
      if( ( basset .eq. 'ablj' ) .or. ( basset .eq. 'bess' ) .or.
     +    ( basset .eq. 'lgsp' ) )then
        if( master )print *, 'sfptab: computing', lsfpr, ' values for',
     +           lastf, ' ', basset, ' functions...'
      else
        call crash( 'SFPTAB', 'Unrecognized basis' )
      end if
c set table of weights for Simpson's rule - table of abscissae not needed
      if( ( basset .eq. 'bess' ) .or. ( basset .eq. 'lgsp' ) )then
        do n = 0, maxn
          sfpwt( n ) = deltak / 3.
          if( ( n .gt. 0 ) .and. ( n .lt. maxn ) )then
            sfpwt( n ) = 2 * sfpwt( n )
            if( mod( n, 2 ) .eq. 1 )sfpwt( n ) = 2 * sfpwt( n )
          end if
        end do
      end if
c select radii to tabulate
      do i = 1, lsfpr
        sfprad( i ) = dble( i - 1 ) / dble( lsfpr - 1 )
        if(
     +    basset .ne. 'ablj' )sfprad( i ) = sfprad( i ) * rgrid( jgrid )
      end do
c compute function values at these radii
      l = 0
      do j = 1, lastf
        n = nsel( j )
        m = msel( j )
        if( basset .eq. 'bess' )k = real( n ) * deltak
        if( basset .eq. 'lgsp' )then
          k = real( n - maxn / 2 ) * deltak
          xcoeff( m, n, 1 ) = nkaln( m, k ) / ( 2 * pi )**2
        end if
        do i = 1, lsfpr
          r = sfprad( i )
          if( basset .eq. 'ablj' )then
            sfplst( 1, l + i ) = ajfun( m, n, r )
            sfplst( 2, l + i ) = ajder( m, n, r )

c            if( mod( i, 50 ) .eq. 0 )print '( 3i5, 3f10.5 )', m, n, i,
c     +  r, sngl( sfplst( 1, l + i ) ), sngl( sfplst( 1, l + i ) )

          else if( basset .eq. 'bess' )then
            x = r * k
            sfplst( 1, l + i ) = bsjfun( m, x )
            sfplst( 2, l + i ) = k * bsjder( m, x )
          else if( ( basset .eq. 'lgsp' ) .and. ( r .gt. 0. ) )then
            x = k * log( r )
            sfplst( 1, l + i ) = cos( x ) / sqrt( r )
            sfplst( 2, l + i ) = sin( x ) / sqrt( r )
          end if
        end do
        l = l + lsfpr
      end do
      if( master )print *, 'sfptab: finished'
      return
      end
