      subroutine frqplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot frequency data from the runXXX.res file
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
c external
      real roundup
c
c local variables
      integer i, ifail, ilast, istrt, j, k, l, nc, ntime
      logical first
      real r, rmn, x, xmax, y, ymax, ymin
c
      call lstrid( 'FRQS', ntime, ilast )
c restart file
      call rewnd
      nc = 0
      istrt = 0
c get first record
      call nextrec( 'FRQS', ifail )
c initialise picture
      call jspage
c set window
      rmn = rinh( jgrid )
      xmax = rgrid( jgrid ) / lscale
c choose first non-zero epicycle frequency as highest
      k = mr + 2
      ymax = roundup( 1.1 * wres( k ) )
      ymax = log10( ymax )
      ymin = ymax - 3
      call jssize( .15, .95, .1, .9 )
      call jscale( 0., xmax, ymin, ymax )
c draw axes
      call jsthik( 3 )
      call jsaxis( 'x', 'Radius', 1 )
      call jslaxs( 'y', 'Frequency', 1 )
      call jsthik( 1 )
c work down file
      do while ( ifail .eq. 0 )
        if( mod( istep, ntime ) .eq. 0 )then
          nc = nc + 1
c work over frequency types
          do j = 1, ma
            first = .true.
            l = j
            if( j .eq. 3 )l = 1
            if( j .eq. 2 )call jsdash( 2, 2, 2, 2 )
            if( j .eq. 3 )call jsdash( 0, 1, 0, 1 )
            do i = l, mr
              r = ( rmn + real( i ) * drfrqs ) / lscale
              k = i + ( j - 1 ) * mr
              y = wres( k )
              if( y .gt. 0. )then
                y = log10( y )
              else
                y = ymin
              end if
              if( first )call jsmove( r, y )
              first = .false.
              call jsline( r, y )
            end do
          end do
          call jsdash( 0, 0, 0, 0 )
        end if
        if( nc .eq. 1 )istrt = istep
        call nextrec( 'FRQS', ifail )
      end do
c label plot
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      if( nc .gt. 1 )then
        call jsbldt( 'time interval' )
        i = ntime
        if( ntime .eq. 1 )i = ( istep - istrt ) / ( nc - 1 )
        x = ts * real( i )
        if( x .lt. 99. )then
          call jsbldf( x, 5, 2 )
        else
          i = x + .5
          call jsbldi( i, 6 )
        end if
      else
        call jsbldt( 't =' )
        x = ts * real( istrt )
        if( x .lt. 99. )then
          call jsbldf( x, 5, 2 )
        else
          i = x + .5
          call jsbldi( i, 6 )
        end if
      end if
      y = 1.05 * ymax - .05 * ymin
      call jsthik( 3 )
      call jswrit( .05 * xmax, y )
      return
      end
