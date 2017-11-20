      subroutine xtplot
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot Toomre's X(r) from the runXXX.res file
c
c Graphics routines are from JSLIB
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      external lxtoom1
      real uofgr
c
c local arrays
      real surf( 150 )
c
c local variables
      integer i, ifail, ilast, istrt, j, k, nc, ntime
      logical first
      real du, r, rmn, srf, u, x, xmax, xmin, y, ymax, ymin
      include 'inc/pi.f'
c
      if( .not. ( p2d .or. p3d ) )return
c choose time interval
      call lstrid( 'DANL', ntime, ilast )
      do i = 1, ncmp
        if( disc( i ) )icmp = i
      end do
c
c initialise picture
      call jspage
c set window
      xmin = 0
      xmax = rtrunc( icmp )
      call yrange( xmin, xmax, ymin, ymax, lxtoom1 )
      call jssize( .15, .95, .1, .9 )
      xmax = rgrid( jgrid ) / lscale
      call jscale( xmin, xmax, ymin, ymax )
c draw axes
      call jsaxis( 'x', 'Radius', 1 )
      call jslaxs( 'y', 'Xm', 1 )
      rmn = rinh( jgrid )
c restart file
      call rewnd
      nc = 0
c
c start of main loop
c
      call nextrec( 'DANL', ifail )
      do while ( ifail .eq. 0 )
        if( nc .eq. 1 )istrt = istep
c density distribution
        if( mod( istep, ntime ) .eq. 0 )then
          nc = nc + 1
          j = 1 - ma
          first = .true.
          do i = 1, mr
            j = j + ma
            surf( i ) = wres( j ) / real( na )
          end do
c frequencies
          call nextrec( 'FRQS', ifail )
          first = .true.
          do i = 2, mr
c epicycle frequency
            k = i + mr
            if( wres( k ) .gt. 0. )then
              r = ( rmn + real( i ) * drfrqs ) / lscale
c get surface density
              u = uofgr( r * lscale )
              j = u
              du = u - real( j )
              srf = surf( j + 1 ) * ( 1. - du ) + surf( j + 2 ) * du
              if( srf .eq. 0. )then
                y = ymin
              else
c compute X
                y = r * wres( k )**2 / ( 2. * pi * srf )
                y = log10( y )
              end if
              if( first )call jsmove( r, y )
              call jsline( r, y )
              first = .false.
            end if
          end do
        end if
        call nextrec( 'DANL', ifail )
        if( istep .gt. ilast )ifail = 1
      end do
c draw theoretical curve
      call jsdash( 2, 2, 2, 2 )
      call jsplt2( lxtoom1 )
      call jsdash( 0, 0, 0, 0 )
c label plot
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      call jsbldt( 'time interval' )
      i = ntime
      if( ( ntime .eq. 1 ) .and.
     +    ( nc .gt. 1 ) )i = ( istep - istrt ) / ( nc - 1 )
      x = ts * real( i )
      if( x .lt. 99. )then
        call jsbldf( x, 5, 2 )
      else
        i = x + .5
        call jsbldi( i, 6 )
      end if
      y = 1.05 * ymax - .05 * ymin
      call jsthik( 3 )
      call jswrit( .95 * xmin + .05 * xmax, y )
      return
      end

      real*8 function lxtoom1( r )
      implicit none
      real*8 r, xtoom
      lxtoom1 = xtoom( r, 1 )
      if( lxtoom1 .gt. 0.d0 )then
        lxtoom1 = log10( lxtoom1 )
      else
        lxtoom1 = 0
      end if
      return
      end
