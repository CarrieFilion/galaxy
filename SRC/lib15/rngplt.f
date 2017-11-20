      subroutine rngplt
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to plot positions of ring particles at requested times
c
c called from analys
c uses JSLIB routines
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
      include 'inc/jscmmn.f'
c
c external
      logical jscren
c
c local arrays
      integer mdim
      parameter ( mdim = 20000 )
      real rad( mdim ), x( mdim ), y( mdim ), z( mdim )
c
c local variables
      character*3 a
      integer i, iend, ifail, iskip, istart, it, j, ka, m, ntime, nw
      logical fl, meriod, positn, skip, tkl
      real rmn, rmx, x1, x2, xl, xmx, y1, y2, yl, ymx, zmx
c
      if( ( jgrid .gt. ngrid ) .or.
     +    ( jgrid .lt. 1 ) )call crash( 'RNGPLT', 'jgrid undefined' )
c
      rmx = rgrid( jgrid ) / lscale
      rmn = rinh( jgrid ) / lscale
      xmx = rmx
      ymx = rmx
      zmx = zm( jgrid ) / lscale
      m = 2
      if( threed )m = 3
      ka = na / nsect
c read in instructions from terminal
      call rewnd
      meriod = .false.
      positn = .false.
      do while ( .not. ( positn .or. meriod ) )
        if( threed )then
          call gtchar(
     +   'Plot projected positions or meriodional, or end?',
     +    a )
        else
          call gtchar( 'Position or end?', a )
        end if
        call jschsz( .3 )
        if( a .eq. 'end' )return
        positn = ( a .eq. 'pos' ) .or. ( a .eq. 'pro' )
        meriod = a .eq. 'mer'
      end do
      call jschsz( .3 )
      call jsthik( 1 )
      call gtintg( 'Enter step number of first frame required', istart )
      iend = 0
      if( .not. jscren( 0 ) )call gtintg(
     +      'Enter step no of last frame required, or 0 for all', iend )
      if( iend .eq. 0 )iend = 1000000
      call gtintg( 'Interval for plotting, or 0 for all?', ntime )
      ntime = max( 1, ntime )
      call gtintg( 'Enter fraction of particles to be plotted', iskip )
      iskip = m * max( 1, iskip )
      tkl = .false.
      if( meriod )then
        x1 = .1
        x2 = x1 + .72
        y1 = .05
        y2 = y1 + .8
        tkl = zmx .gt. 1.
        xl = x2 + .01
        yl = .5 * ( y1 + y2 )
      end if
c work through file
      call jschsz( .2 )
      skip = .true.
      do while ( skip )
        call nextrec( 'RNGS', ifail )
        if( ( ifail .ne. 0 ) .or. ( istep .gt. iend ) )then
          call jschsz( .3 )
          return
        end if
        skip = ( mod( istep, ntime ) .ne. 0 ) .or. ( istep .lt. istart )
        if( .not. skip )then
          nw = m * mr * ma
          if( nw / iskip .gt. mdim )call space(
     +                       mdim, nw / iskip, 'Work arrays', 'RNGPLT' )
c write heading
          call jspage
          call jssize( 0., 1., 0., 1. )
          call jscale( 0., 1., 0., 1. )
          if( ipic .eq. 1 )then
            call jsbldt( 'Run no' )
            call jsbldi( irun, 6 )
            call jswrit( .3, .95 )
          end if
          call jsbldt( 't=' )
          it = time + .5
          fl = abs( time - real( it ) ) .gt. .01
          if( fl )then
            call jsbldf( time, 5, 1 )
          else
            call jsbldi( it, 3 )
          end if
          call jswrit( .75, .83 )
c set up frame
          call jssize( .05, .95, 0., .9 )
          call jsescl( -rmx, rmx, -rmx, rmx )
          if( rmn .gt. 0. )call jscirc( 0., 0., rmn )
          if( p2d .or. p3d )call jscirc( 0., 0., rmx )
c calculate coordinates
          j = 0
          do i = 1, nw, iskip
            if( wres( i ) .ne. 0. )then
              j = j + 1
              x( j ) = wres( i ) / lscale
              y( j ) = wres( i + 1 ) / lscale
              z( j ) = wres( i + 2 ) / lscale
              rad( j ) = sqrt( x( j )**2 + y( j )**2 )
            end if
          end do
c plot points
          if( .not. jscren(0) )print '( '' Plotting run no'', i5, ' //
     +      ' '' step no'', 2i6, '' particles'' )', irun, istep, j
          if( positn )then
            call posplt( x, y, z, j, .true., rmx, zmx )
c meriodional projection
          else if( meriod )then
            call jssize( x1, x2, y1, y2 )
            call jsescl( 0., xmx, -zmx, zmx )
            call jsaxis( 'x', 'r', 1 )
            if( tkl )then
              call jsaxis( 'y', 'z', 1 )
            else
              call jsaxis( 'y', 'z', 0 )
            end if
            call jsymbl( rad, z, j, 0 )
          end if
        end if
        skip = .true.
      end do
      end
