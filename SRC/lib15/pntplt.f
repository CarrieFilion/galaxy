      subroutine pntplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Plots positions or velocities of the particles in the runXXX.res file
c   with suitable annotation
c
c Called from ANALYS
c Graphics routines are from JSPLOT
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      include 'inc/jscmmn.f'
c
      include 'inc/model.f'
c
c externals
      logical gtlogl
      real grofu, roundup
c
c local arrays
      integer ibuf
      parameter ( ibuf = 500000 )
      real, allocatable :: rad(:), u(:), v(:), vw(:), x(:), y(:), z(:)
c
c local variables
      character*3 a
      integer i, iend, ifail, ip, iskip, istart, istrd, it, j, ka, sp
      integer nextc, ntime
      logical all, fl, meriod, need, positn, recen, tkl, veloc
      real r, rmn, rmx, xl, xmx, x1, x2, x3, x4, x5, x6, yl, ymx, yu, yv
      real yw, y1, y2, zmx, smx, xc, yc, zc
c
      ka = na / nsect
c read in instructions from terminal
      call rewnd
      call jsthik( 1 )
      meriod = .false.
      positn = .false.
      veloc = .false.
      do while ( .not. ( positn .or. meriod .or. veloc ) )
        if( threed )then
          call gtchar(
     +   'Plot projected positions, meriodional or velocities, or end?',
     +    a )
        else
          call gtchar( 'Position or velocity plot, or end?', a )
        end if
        call jschsz( .3 )
        if( a .eq. 'end' )return
        positn = ( a .eq. 'pos' ) .or. ( a .eq. 'pro' )
        meriod = a .eq. 'mer'
        veloc = a .eq. 'vel'
      end do
c
      r = 0
      ip = 0
      do i = 1, ncmp
        if( nsp( i ) .gt. 0 )then
          ip = ip + 1
c identify grid with largest rmax
          j = igrd( i )
          if( rgrid( j ) .gt. r )then
            ka = j
            r = rgrid( j )
          end if
        end if
      end do
c switch to grid with largest rgrid
      call switch( ka )
      if( rngs )ip = ip - 1
      if( ip .le. 0 )call crash( 'PNTPLT', 'No active pop selected' )
c select which pop to plot, or all
      if( ip .gt. 1 )then
        all =
     +     gtlogl( 'Multiple active components - do you want them all' )
        if( .not. all )then
          call selpop( sp )
          call switch( igrd( sp ) )
          print *, 'population', sp, ' selected'
        end if
      else
        all = .true.
        sp = 1
      end if
c
      rmx = 0
      if( twod )then
        rmx = rgrid( jgrid ) / lscale
c        print *, 'The grid edge is at', rmx
c        call gtchar(
c     +        'Do you want to draw the circle at another radius? (y/n)',
c     +        a( 1:1 ) )
c        if( a( 1:1 ) .eq. 'y' )then
c          call gtreal( 'Enter new radius', rmn )
c          rmx = min( rmx, rmn )
c        end if
        if( uoffs .eq. 0. )then
          rmn = 0
        else
          rmn = rinh( jgrid ) / lscale
        end if
        xmx = rmx
        ymx = rmx
      else
c other grid types
        rmn = 0
        rmx = rgrid( jgrid ) / lscale
c$$$        if( ( .not. all ) .and. ( rmx .gt. 2. * rtrunc( sp ) ) )then
c$$$          rmx = 1.5 * rtrunc( sp )
c$$$        else
          print *, 'rmx =', rmx
          if( gtlogl( 'Do you want to change rmx' ) )call gtreal(
     +                                 'Enter rmx in disk scales', rmx )
c$$$        end if
        xmx = rmx
        ymx = rmx
        zmx = zm( jgrid ) / lscale
        zmx = min( zmx, rmx )
        if( positn )then
          print *, 'zmx =', zmx
          if( gtlogl( 'Do you want to change zmx' ) )call gtreal(
     +                                 'Enter zmx in disk scales', zmx )
        end if
        if( p3a )rmn = grofu( 0. ) / lscale
        if( p3d )rmx = max( xmx, zmx )
        if( c3d )then
          xmx = xm / lscale
          ymx = ym / lscale
          rmx = max( xmx, ymx, zmx )
        end if
      end if
c
      call lstrid( 'PNTS', ntime, iend )
      call gtreal( 'Enter time of first frame required', r )
      istart = nint( r / ts )
      call gtintg( 'Enter fraction of particles to be plotted', iskip )
      iskip = ncoor * max( 1, iskip )
      tkl = .false.
c up set frame boundaries
      if( positn )then
c position space
        if( twod )then
          if( c2d )then
            xl = .2
            yl = .90
          else
            xl = .75
            yl = .83
          end if
        else
          smx = max( xmx + zmx, ymx + zmx, 1.2 * xmx )
          smx = .85 / smx
          tkl = zmx .gt. 1.
          xl = .11 + smx * xmx
          yl = .1 + smx * ( ymx + .5 * zmx )
        end if
      else if( meriod )then
        x1 = .1
        x2 = x1 + .72
        y1 = .15
        y2 = y1 + .8
        tkl = zmx .gt. 1.
        xl = x2 + .01
        yl = .5 * ( y1 + y2 )
      else if( veloc )then
c phase space
        if( twod )then
          x1 = .05
          x2 = x1 + .4
          x3 = x2 + .05
          x4 = x3 + .4
        else
          x1 = .05
          x2 = x1 + .28
          x3 = x2 + .05
          x4 = x3 + .28
          x5 = x4 + .05
          x6 = x5 + .28
        end if
        xl = x4 - .1
        yl = .83
      end if
c determine whether to shift the frame with the grid
      if( centrd )then
        recen = gtlogl( 'Move frame with the grid' )
      else
        recen = .false.
      end if
      if( recen )then
        call nextrec( 'XCEN', ifail )
        xc = wres( 1 )
        yc = wres( 2 )
        zc = wres( 3 )
      else
        xc = 0
        yc = 0
        zc = 0
      end if
c allocate space
      allocate ( rad( ibuf ) )
      allocate ( x( ibuf ) )
      allocate ( y( ibuf ) )
      if( threed )allocate ( z( ibuf ) )
      if( veloc )then
        allocate ( u( ibuf ) )
        allocate ( v( ibuf ) )
        if( threed )allocate ( vw( ibuf ) )
      end if
c work through file
      call jschsz( .2 )
      gvfac = 1. / ( lscale * ts )
      call rewnd
      call nextrec( 'PNTS', ifail )
      nextc = 1
      do while ( ifail .eq. 0 )
        if( ( mod( istep, ntime ) .eq. 0 ) .and.
     +      ( istep .ge. istart ) )then
c calculate coordinates
          if( veloc )then
            yu = .1
            yv = .1
            yw = .1
          end if
          j = 0
          do i = 1, mr, iskip
            if( j .ge. ibuf )call crash( 'PNTPLT',
     +                                         'work arrays too small' )
            if( twod )then
              if( wres( i ) .ne. 0. )then
                j = j + 1
                if( positn )then
                  x( j ) = wres( i ) / lscale
                  y( j ) = wres( i + 1 ) / lscale
                  if( x( j )**2 + y( j )**2 .gt. rmx**2 )j = j - 1
                else
                  r = sqrt( wres( i )**2 + wres( i + 1 )**2 )
                  u( j ) = ( wres( i ) * wres( i + 2 ) +
     +                       wres( i + 1 ) * wres( i + 3 ) ) * gvfac / r
                  v( j ) = ( wres( i ) * wres( i + 3 ) -
     +                       wres( i + 1 ) * wres( i + 2 ) ) * gvfac / r
                  rad( j ) = r / lscale
                  yu = max( yu, abs( u( j ) ) )
                  yv = max( yv, abs( v( j ) ) )
                end if
              end if
            else
c allow for coordinates in either order
              j = j + 1
              x( j ) = wres( i ) - xc
              y( j ) = wres( i + 1 ) - yc
              z( j ) = wres( i + 2 ) - zc
              if( positn )then
                x( j ) = x( j ) / lscale
                y( j ) = y( j ) / lscale
                rad( j ) = sqrt( x( j )**2 + y( j )**2 )
              else if( meriod )then
                rad( j ) = sqrt( x( j )**2 + y( j )**2 ) / lscale
              else if( veloc )then
                rad( j ) = sqrt( x( j )**2 + y( j )**2 )
                u( j ) = ( x( j ) * wres( i + 3 ) +
     +                     y( j ) * wres( i + 4 ) ) * gvfac / rad( j )
                v( j ) = ( x( j ) * wres( i + 4 ) -
     +                     y( j ) * wres( i + 3 ) ) * gvfac / rad( j )
                vw( j ) = wres( i + 5 ) * gvfac
                rad( j ) = rad( j ) / lscale
                yu = max( yu, abs( u( j ) ) )
                yv = max( yv, abs( v( j ) ) )
                yw = max( yw, abs( vw( j ) ) )
              end if
c pick up and remove particle population flag
              ip = nint( z( j ) / zshift )
              z( j ) = ( z( j ) - real( ip ) * zshift ) / lscale
              ip = ip + 1
c skip unwanted particles, not selected pop
              need = all .or. ( ip .eq. sp )
c radius too large
              if( positn )then
                need = need .and. ( x( j )**2 + y( j )**2 .le. rmx**2 )
              else
                need = need .and. ( rad( j ) .le. rmx )
              end if
c z height too large
              if( threed )then
                need = need .and. ( abs( z( j ) ) .le. zmx )
              end if
              if( .not. need )j = j - 1
            end if
          end do
          if( j .gt. 0 )then
            print '( ''Plotting run'', i5, '' time'', f8.1, i8, ' //
     +            ' '' particles'' )', irun, time, j
c write heading
            call jspage
            call jssize( 0., 1., 0., 1. )
            call jscale( 0., 1., 0., 1. )
            if( ipic .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 6 )
              call jswrit( xl, .97 )
            end if
            call jsbldt( 't=' )
            it = time + .5
            fl = abs( time - real( it ) ) .gt. .01
            if( fl )then
              call jsbldf( time, 5, 1 )
            else
              if( it .gt. 9999 )then
                call jsbldi( it, 5 )
              else if( it .gt. 999 )then
                call jsbldi( it, 4 )
              else
                call jsbldi( it, 3 )
              end if
            end if
            call jswrit( xl, yl )
c position space
            if( positn )then
              call posplt( x, y, z, j, .true., rmx, zmx )
            else if( meriod )then
c meriodional projection
              call jssize( x1, x2, y1, y2 )
              if( scale_set )then
                call jsescl( 0., xmx * unit_L,
     +                       -zmx * unit_L, zmx * unit_L )
                call jsaxis( 'x', 'R (kpc)', 1 )
                if( tkl )then
                  call jsaxis( 'y', 'z (kpc)', 1 )
                else
                  call jsaxis( 'y', 'z (kpc)', 0 )
                end if
                call jsescl( 0., xmx, -zmx, zmx )
              else
                call jsescl( 0., xmx, -zmx, zmx )
                call jsaxis( 'x', 'R', 1 )
                if( tkl )then
                  call jsaxis( 'y', 'z', 1 )
                else
                  call jsaxis( 'y', 'z', 0 )
                end if
              end if
              call jsymbl( rad, z, j, 0 )
            else if( veloc )then
c phase space
              yu = roundup( yu )
              yv = roundup( yv )
              yw = roundup( yw )
c force equal scaling and reasonable range
              yu = max( yu, yv, yw )
              yu = min( yu, 2. )
              yv = yu
              yw = yu
c plot velocity components
              call jssize( x1, x2, .1, .9 )
              if( scale_set )then
                call jscale( rmn * unit_L, rmx * unit_L,
     +                       -yu * unit_V, yu * unit_V )
                call jsaxis( 'x', 'R (kpc)', 1 )
                call jsaxis( 'y', 'u (km/s)', 1 )
                call jscale( rmn, rmx, -yu, yu )
              else
                call jscale( rmn, rmx, -yu, yu )
                call jsaxis( 'x', 'R', 1 )
                call jsaxis( 'y', 'u', 1 )
              end if
              call jsymbl( rad, u, j, 0 )
c
              call jssize( x3, x4, .1, .9 )
              if( scale_set )then
                call jscale( rmn * unit_L, rmx * unit_L,
     +                       -yu * unit_V, yv * unit_V )
                call jsaxis( 'x', 'R (kpc)', 1 )
                call jsaxis( 'y', 'v (km/s)', 1 )
                call jscale( rmn, rmx, -yu, yv )
              else
                call jscale( rmn, rmx, -yu, yv )
                call jsaxis( 'x', 'r', 1 )
                call jsaxis( 'y', 'v', 1 )
              end if
              call jsymbl( rad, v, j, 0 )
c
              if( threed )then
                call jssize( x5, x6, .1, .9 )
                if( scale_set )then
                  call jscale( rmn * unit_L, rmx * unit_L,
     +                         -yw * unit_V, yw * unit_V )
                  call jsaxis( 'x', 'R (kpc)', 1 )
                  call jsaxis( 'y', 'w (km/s)', 1 )
                  call jscale( rmn, rmx, -yw, yw )
                else
                  call jscale( rmn, rmx, -yw, yw )
                  call jsaxis( 'x', 'R', 1 )
                  call jsaxis( 'y', 'w', 1 )
                end if
                call jsymbl( rad, vw, j, 0 )
              end if
            end if
          end if
        end if
        istrd = istep
c get next center if needed
        if( recen )then
          do while ( istep .lt. nextc )
            call nextrec( 'XCEN', ifail )
            if( ifail .ne. 0 )istep = nextc
          end do
          xc = wres( 1 )
          yc = wres( 2 )
          zc = wres( 3 )
        end if
c get next data record
        call nextrec( 'PNTS', ifail )
        istrd = istep - istrd
        nextc = istep + istrd
        if( istep .gt. iend )ifail = 1
      end do
      call jsthik( 3 )
      call jschsz( .3 )
      return
      end
