      subroutine satanl
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c part of the analysis software
c it provides a number of interactive options for displaying the
c   the orbits of the host and satellite
c   OR the motions of the grid center(s) when they are allowed to move
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/comprns.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      real roundup
c
c local arrays
      integer, allocatable :: is( : )
      real, allocatable :: sat( :, : ), xs( : ), ys( : ), zs( : )
      real, allocatable :: xc( :, : ), vc( :, : )
c
c local variables
      character*4 type
      integer i, iamx, ifail, iopt, isub, j, k, kgrid, krun, m, me
      integer nsteps
      logical tkl, twog
      real actmas, amx, d, fac, ri, smx, t, tlast, xl, x1, x2, x3, x4
      real xk( 2 ), xmin, xmax, y, yl, y1, y2, y3, y4, yk( 2 ), ylast
      real ymax, ymin, yrange, zmin, zmax
      equivalence ( i, ri )
c
c select type of data
      type = ' '
      do while ( ( type .ne. 'SATE' ) .and. ( type .ne. 'XCEN' ) )
        call gtchar( 'SATE or XCEN data', type )
        call uppercase( type )
      end do
      sate = type .eq. 'SATE'
c check whether requested data exists
      call rewnd
      call nextrec( type, ifail )
      if( ifail .ne. 0 )then
        print *, 'No data of requested type found'
        go to 6
      end if
      me = mr
      if( .not. sate )me = mr * ma
c count number of records available
      nsteps = 0
      do while( ifail .eq. 0 )
        nsteps = nsteps + 1
        call nextrec( type, ifail )
      end do
c allocate space
      allocate ( sat( me, nsteps ) )
      allocate ( xs( nsteps ), ys( nsteps ), zs( nsteps ) )
      allocate ( is( nsteps ), xc( 3, nsteps ), vc( 3, nsteps ) )
c read in data
      krun = 1
    5 call rewnd
      j = 0
      ifail = 0
      do while( ifail .eq. 0 )
        j = j + 1
c position and drift of CoM
        if( sate )then
          call nextrec( 'INTG', ifail )
          if( ifail .eq. 0 )then
            do i = 1, 3
              xc( i, j ) = 0
              vc( i, j ) = 0
            end do
            m = 5
            actmas = 0
            do k = 1, ncmp
              if( cmpmas( k ) .gt. 0. )then
                if( uqmass )then
                  fac = wres( m )
                else
                  ri = wres( m )
                  fac = i
                end if
                actmas = actmas + fac
                do i = 1, 3
                  xc( i, j ) = xc( i, j ) + wres( i + m ) * fac
                  vc( i, j ) = vc( i, j ) + wres( i + m + 3 )
                end do
              end if
              m = m + 15
            end do
            do i = 1, 3
              xc( i, j ) = xc( i, j ) / actmas
            end do
          end if
        end if
c get next main record
        call nextrec( type, ifail )
        if( j .gt. nsteps )ifail = 1
        if( ifail .eq. 0 )then
          is( j ) = istep
          if( sate )then
            do i = 1, me
              sat( i, j ) = wres( i )
            end do
            actmas = 0
            do k = 1, ncmp
              actmas = actmas + cmpmas( k )
            end do
            fac = actmas + mptbr
            do i = 1, 3
              xc( i, j ) = ( actmas * xc( i, j ) +
     +                       mptbr * wres( i ) ) / fac
              vc( i, j ) = vc( i, j ) + mptbr * wres( i + 3 )
            end do
          else
            do i = 1, me
              sat( i, j ) = wres( i ) / lscale
            end do
            do i = 1, 3
              xc( i, j ) = 0
              vc( i, j ) = 0
            end do
          end if
          if( exsphr )sat( 10, j ) = wres( 10 ) +
     +       .5 * mptbr * ( wres( 4 )**2 + wres( 5 )**2 + wres( 6 )**2 )
        end if
      end do
      if( krun .eq. 1 )tlast = is( nsteps ) * ts
      if( extbar )then
c smooth using running average
        do j = 1, nsteps
          xs( j ) = sat( 3, j )
        end do
        k = 50
c        call gtintg( 'Enter k', k )
        k = min( k, nsteps )
        call bcsmth( xs, nsteps, k )
        do j = 1, nsteps
          sat( 3, j ) = xs( j )
        end do
cc compute acceleration from derivative of pattern speed
c        do j = 1, nsteps
c          xs( j ) = ts * real( is( j ) )
c          ys( j ) = sat( 2, j )
c        end do
c        k = 50
cc        call gtintg( 'Enter k', k )
c        k = min( k, nsteps )
c        do j = 1, nsteps
c          i = max( 1, j - k / 2 )
c          i = min( i, nsteps - k + 1 )
c          call linlsq( xs( i ), ys( i ), k, d, y1, y2, y3 )
c          sat( 3, j ) = d
c        end do
      end if
c rescale and find largest
      if( extbar )then
        amx = 0
        do j = 1, nsteps
          sat( 3, j ) = sat( 3, j ) / mptbr
          if( -sat( 3, j ) .gt. amx )then
            amx = -sat( 3, j )
            iamx = j
          end if
        end do
        print *, 'abs max accn =', amx, ' at moment', iamx
        print *, 'when omega_b =', sat( 2, iamx )
      end if
      if( krun .eq. 1 )then
        print *, nsteps, ' records read in to time', tlast
        call jssize( .2, .95, .15, .95 )
        xmax = roundup( .99 * tlast)
      end if
      if( sate )then
c choose option for satellite data
    1   if( krun .eq. 1 )then
          call gtintg( 'Enter option, -1 to end or 0 for help', iopt )
          if( iopt .lt. 0 )go to 6
          xmax = roundup( .99 * tlast)
          if( iopt .eq. 0 )then
            if( exsphr .or. exgnrc )then
              print *, '1 - z position of satellite'
              print *, '2 - z velocity of satellite'
              print *, '3 - acc of satellite'
              print *, '5 - total E of satellite'
              print *, '6 - orbit of satellite'
              print *, '7 - R(t)'
              print *, '8 - L_z(t)'
              if( m .gt. 10 )print *, '9 - M(t)'
            else if( extbar )then
              print *, '1 - phase of bar'
              print *, '2 - pattern speed of bar'
              print *, '3 - accn of bar'
              print *, '4 - speed vs accn of bar'
            else
              call crash( 'SATANL', 'Unrecognized perturber' )
            end if
            go to 1
          end if
          if( iopt .gt. 9 )go to 1
          if( exsphr .and. ( iopt .eq. 4 ) )go to 1
          if( ( iopt .eq. 9 ) .and. ( m .lt. 11 ) )go to 1
          if( extbar .and. ( iopt .gt. 4 ) )go to 1
        end if
c projected orbit plots
        if( iopt .eq. 6 )then
          do i = 1, nsteps
            xs( i ) = sat( 1, i ) - xc( 1, i )
            ys( i ) = sat( 2, i ) - xc( 2, i )
            zs( i ) = sat( 3, i ) - xc( 3, i )
          end do
          if( krun .eq. 1 )then
            xmax = 0
            ymax = 0
            zmax = 0
            do i = 1, nsteps
              xmax = max( xmax, abs( xs( i ) ) )
              ymax = max( ymax, abs( ys( i ) ) )
              zmax = max( zmax, abs( zs( i ) ) )
            end do
            xmax = max( xmax, ymax )
            xmax = roundup( xmax )
            ymax = xmax
            zmax = roundup( zmax )
c set frame markers
            smx = max( xmax + zmax, ymax + zmax, 1.2 * xmax )
            smx = .85 / smx
            x1 = .1
            x2 = x1 + smx * xmax
            x3 = x2
            x4 = x3 + smx * zmax
            y1 = .1
            y2 = y1 + smx * ymax
            y3 = y2
            y4 = y3 + smx * zmax
            tkl = zmax .gt. 1.
            xl = .11 + smx * xmax
            yl = .1 + smx * ( ymax + .5 * zmax )
            call jspage
          end if
          call jssize( 0., 1., 0., 1. )
          call jsescl( 0., 1., 0., 1. )
          if( nruns .eq. 1 )then
            call jsbldt( 'Run no' )
            call jsbldi( irun, 6 )
            call jswrit( xl, .97 )
          else
            xk( 1 ) = .6
            xk( 2 ) = .65
            yk( 1 ) = .5 + .1 * real( krun )
            yk( 2 ) = yk( 1 )
            call jscoin( xk, yk, 2, krun )
            call jsbldi( irun, 4 )
            call jswrit( .66, yk( 1 ) - .02 )
          end if
c x-y projection
          call jssize( x1, x2, y1, y2 )
          call jsescl( -xmax, xmax, -ymax, ymax )
          if( krun .eq. 1 )then
            call jsaxis( 'x', 'x', 1 )
            call jsaxis( 'y', 'y', 1 )
          end if
          call jscoin( xs, ys, nsteps, krun )
          if( zmax .gt. 0.01 * min( xmax, ymax ) )then
c full 3-D plot
            call jssize( x3, x4, y1, y2 )
            call jsescl( -zmax, zmax, -ymax, ymax )
            if( krun .eq. 1 )then
              if( tkl )then
                call jsaxis( 'x', 'z', 1 )
              else
                call jsaxis( 'x', 'z', 0 )
              end if
              call jsaxis( 'y', ' ', 0 )
            end if
            call jscoin( zs, ys, nsteps, krun )
c
            call jssize( x1, x2, y3, y4 )
            call jsescl( -xmax, xmax, -zmax, zmax )
            if( krun .eq. 1 )then
              call jsaxis( 'x', ' ', 0 )
              if( tkl )then
                call jsaxis( 'y', 'z', 1 )
              else
                call jsaxis( 'y', 'z', 0 )
              end if
            end if
            call jscoin( xs, zs, nsteps, krun )
          end if
          call jssize( .2, .95, .15, .95 )
          xmax = roundup( .99 * tlast)
          go to 4
c radius in the (x,y) plane
        else if( iopt .eq. 7 )then
          do i = 1, nsteps
            xs( i ) = ts * real( is( i ) )
            ys( i ) = sqrt( ( sat( 1, i ) - xc( 1, i ) )**2 +
     +                      ( sat( 2, i ) - xc( 2, i ) )**2 )
          end do
          if( krun .eq. 1 )then
            ymax = 0
            ymin = 1.e10
            do i = 1, nsteps
              ymin = min( ymin, abs( ys( i ) ) )
              ymax = max( ymax, abs( ys( i ) ) )
            end do
            i = roundup( ymax - ymin )
            i = max( i, 2 )
            ymin = int( ymin - .5 )
            ymax = max( ymin + i, ymax )
            call jspage
            call jscale( 0., xs( nsteps ), ymin, ymax )
            call jsaxis( 'x', 't', 1 )
            call jsaxis( 'y', 'R', 1 )
          end if
          call jscoin( xs, ys, nsteps, krun )
          go to 4
        end if
c convert option number to array subscript
        if( exsphr )then
          if( iopt .eq. 1 )isub = 3
          if( iopt .eq. 2 )isub = 6
          if( iopt .eq. 3 )isub = 9
          if( iopt .eq. 5 )isub = 10
          if( iopt .eq. 9 )isub = 11
        else
          isub = iopt
          if( extbar .and. ( iopt .eq. 4 ) )isub = 3
        end if
c set scale
        if( krun .eq. 1 )then
          ymax = 0
          ymin = 0
          do i = 1, nsteps
            if( iopt .ne. 8 )then
              y = sat( isub, i )
            else
              y = sat( 1, i ) * sat( 5, i ) - sat( 2, i ) * sat( 4, i )
            end if
            ymax = max( ymax, y )
            ymin = min( ymin, y )
          end do
          ymax = roundup( ymax )
          ymin = -roundup( -ymin )
          yrange = ymax - ymin
          if( yrange .eq. 0. )go to 1
          call jspage
          if( exsphr )then
            call jscale( 0., xmax, ymin, ymax )
            call jsaxis( 'x', 'time', 1 )
            if( iopt .eq. 1 )call jsaxis( 'y', 'z_{sat}-z_{disk}', 1 )
            if( iopt .eq. 2 )call jsaxis( 'y', 'v_{sat}', 1 )
            if( iopt .eq. 3 )call jsaxis( 'y', 'a_{sat}', 1 )
            if( iopt .eq. 4 )call jsaxis( 'y', 'zshift', 1 )
            if( iopt .eq. 5 )call jsaxis( 'y', '(Total E)_{sat}', 1 )
            if( iopt .eq. 8 )call jsaxis( 'y', 'L_{z}(sat)', 1 )
            if( iopt .eq. 9 )call jsaxis( 'y', 'M_{central}', 1 )
          else if( extbar )then
            if( iopt .eq. 4 )then
c              xmax = sat( 2, 1 )
              xmax = 2
              call jscale( 0., xmax, -1.2, 0. )
              call jsaxis( 'x', 'bar pattern speed', 1 )
              call jsaxis( 'y', 'bar acceleration/Mbar', 1 )
c              call jscale( 0., xmax, -.4 * mptbr, 0. )
            else
              call jscale( 0., xmax, ymin, ymax )
              call jsaxis( 'x', 'time', 1 )
            end if
            if( iopt .eq. 1 )call jsaxis( 'y', 'bar phase', 1 )
            if( iopt .eq. 2 )call jsaxis( 'y', 'bar pattern speed', 1 )
            if( iopt .eq. 3 )call jsaxis( 'y', 'bar acceleration', 1 )
          end if
        end if
c plot result
        if( nruns .eq. 1 )then
          call jsbldt( 'Run' )
          call jsbldi( irun, 4 )
          call jswrit( .75 * xmax, .9 * ymax + .1 * ymin )
        end if
        if( extbar .and. ( iopt .eq. 4 ) )then
          do i = 1, nsteps
c          call jsymbl( sat( 2, i ), sat( 3, i ), 1, 0 )
            if( i .eq. 1 )call jsmove( sat( 2, i ), sat( 3, i ) )
            call jsline( sat( 2, i ), sat( 3, i ) )
          end do
        else
          if( iopt .eq. 8 )then
            ylast = sat( 1, 1 ) * sat( 5, 1 ) - sat( 2, 1 ) *sat( 4, 1 )
          else
            ylast = sat( isub, 1 )
          end if
          tlast = 0
          call jsmove( 0., ylast )
          do i = 1, nsteps
            t = ts * real( is( i ) )
            if( iopt .eq. 8 )then
              y =
     +       ( sat( 1, i ) - xc( 1, i ) ) * ( sat( 5, i ) - vc( 2, i ) )
     +     - ( sat( 2, i ) - xc( 2, i ) ) * ( sat( 4, i ) - vc( 1, i ) )
            else
              y = sat( isub, i )
            end if
            if( ( abs( y - ylast ) .gt. .005 * yrange ) .or.
     +           ( t - tlast .gt. 1. ) )then
              call jsline( t, y )
              tlast = t
              ylast = y
            end if
          end do
          call jsline( t, y )
        end if
        go to 4
c
      else
c choose option for xcen data
    2   if( krun .eq. 1 )then
          call gtintg( 'Enter option, -1 to end or 0 for help', iopt )
          if( iopt .lt. 0 )go to 6
          if( iopt .eq. 0 )then
            print *, '1 - orbits of centres'
            print *, '2 - distance between centres'
            print *, '3 - radial displacement of centres'
            go to 2
          end if
          if( iopt .gt. 3 )go to 4
          kgrid = ngrid
          if( pertbn .and. ( exsphr .or. exgnrc ) )kgrid = kgrid + 1
          if( twogrd )then
            call switch( 2 )
            twog = ( .not. hybrid ) .and. ( .not. dr3d )
            if( .not. twog )kgrid = 1
          end if
        end if
c
        if( iopt .eq. 1 )then
          if( krun .eq. 1 )then
            xmin = 100
            xmax = -xmin
            ymin = 100
            ymax = -xmin
            zmin = 100
            zmax = -xmin
            do i = 1, nsteps
              k = 0
              do j = 1, kgrid
                xmin = min( xmin, sat( k + 1, i ) )
                xmax = max( xmax, sat( k + 1, i ) )
                ymin = min( ymin, sat( k + 2, i ) )
                ymax = max( ymax, sat( k + 2, i ) )
                if( threed )then
                  zmin = min( zmin, sat( k + 3, i ) )
                  zmax = max( zmax, sat( k + 3, i ) )
                end if
                k = k + ndimen
              end do
            end do
            xmax = 1.05 * max( xmax, -xmin )
            ymax = max( ymax, -ymin )
            zmax = max( zmax, -zmin )
c force square frames
            xmax = max( xmax, ymax, zmax )
            ymax = xmax
            zmax = xmax
c draw and mark axes
            if( twod )then
              x1 = .1
              x2 = x1 + .8
              y1 = .1
              y2 = y1 + .8
              xl = -.7 * xmax
              yl = 1.05 * ymax
            else
              x1 = .1
              x2 = x1 + .4
              x3 = x2
              x4 = x3 + .4
              y1 = .1
              y2 = y1 + .4
              y3 = y2
              y4 = y3 + .4
              xl = .15 * x2
              yl = 2.05 * y2
            end if
            call jspage
c mark header
            if( twog )then
              call jsbldt( 'Grid centres' )
            else
              call jsbldt( 'Grid centre' )
            end if
            call jsbldt( '(model units)' )
          end if
          call jssize( 0., 1., 0., 1. )
          call jscale( 0., 1., 0., 1. )
          if( nruns .eq. 1 )then
            call jsbldt( 'run no' )
            call jsbldi( irun, 4 )
            call jswrit( xl, yl )
          else
            if( krun .eq. 1 )call jswrit( .6, .85 )
            xk( 1 ) = .6
            xk( 2 ) = .65
            yk( 1 ) = .6 + .05 * real( krun )
            yk( 2 ) = yk( 1 )
            call jscoin( xk, yk, 2, krun )
            call jsbldi( irun, 4 )
            call jswrit( .66, yk( 1 ) - .01 )
          end if
c plot data
          k = 0
          do j = 1, kgrid
            call jssize( x1, x2, y1, y2 )
            call jsescl( -xmax, xmax, -ymax, ymax )
            if( j .eq. 1 .and. krun .eq. 1 )then
              call jsaxis( 'x', 'x', 1 )
              call jsaxis( 'y', 'y', 1 )
            end if
            do i = 1, nsteps
              xs( i ) = sat( k + 1, i )
              ys( i ) = sat( k + 2, i )
              if( threed )zs( i ) = sat( k + 3, i )
            end do
            call jscoin( xs, ys, nsteps, krun )
            k = k + ndimen
            if( threed )then
              call jssize( x3, x4, y1, y2 )
              call jsescl( -zmax, zmax, -ymax, ymax )
              if( j .eq. 1 .and. krun .eq. 1 )then
                call jsaxis( 'x', 'z', 1 )
                call jsaxis( 'y', ' ', 0 )
              end if
              call jscoin( zs, ys, nsteps, krun )
              call jssize( x1, x2, y3, y4 )
              call jsescl( -xmax, xmax, -zmax, zmax )
              if( j .eq. 1 .and. krun .eq. 1 )then
                call jsaxis( 'x', ' ', 0 )
                call jsaxis( 'y', 'z', 1 )
              end if
              call jscoin( xs, zs, nsteps, krun )
            end if
          end do
c
        else if( iopt .eq. 2 )then
c
          if( krun .eq. 1 )then
            if( kgrid .ne. 2 )then
              print *, 'Option requires 2 independently moving grids'
              go to 2
            end if
            ymax = -xmin
          end if
          do i = 1, nsteps
            d = 0
            do j = 1, ndimen
              d = d + ( sat( j, i ) - sat( j + ndimen, i ) )**2
            end do
            d = sqrt( d )
            if( krun .eq. 1 )ymax = max( ymax, d )
            xs( i ) = ts * real( is( i ) )
            ys( i ) = d
          end do
          if( krun .eq. 1 )then
            ymax = roundup( ymax )
            xmax = xs( nsteps )
            call jspage
c plot data
            call jssize( .15, .95, .15, .95 )
            call jscale( 0., xmax, 0., ymax )
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', 'distance', 1 )
          end if
          call jscoin( xs, ys, nsteps, krun )
c mark header
          if( nruns .eq. 1 )then
            call jsbldt( 'Distance between grid centres' )
            call jsbldt( 'Run no' )
            call jsbldi( irun, 4 )
            call jswrit( .1 * xmax, 1.1 * ymax )
          else
            xk( 1 ) = .7 * xmax
            xk( 2 ) = .75 * xmax
            yk( 1 ) = ( .95 - .05 * real( krun ) ) * ymax
            yk( 2 ) = yk( 1 )
            call jscoin( xk, yk, 2, krun )
            call jsbldi( irun, 4 )
            call jswrit( xk( 2 ) + .01 * xmax, yk( 1 ) - .01 * ymax )
          end if
c
        else if( iopt .eq. 3 )then
c
c find maximum radius
          if( krun .eq. 1 )then
            ymax = 0
            do k = 1, kgrid
              do i = 1, nsteps
                d = 0
                do j = 1, ndimen
                  d = d + sat( j + ( k - 1 ) * ndimen, i )**2
                end do
                ymax = max( ymax, d )
              end do
            end do
c set up frame
            ymax = roundup( sqrt( ymax ) )
            xmax = ts * real( is( nsteps ) )
            call jspage
            call jssize( .15, .95, .15, .95 )
            call jscale( 0., xmax, 0., ymax )
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', 'distance', 1 )
c mark header
            if( twog )then
              call jsbldt( 'Radial displacement of grid centres' )
            else
              call jsbldt( 'Radial displacement of grid centre' )
            end if
          end if
          if( nruns .eq. 1 )then
            call jsbldt( 'Run no' )
            call jsbldi( irun, 4 )
            call jswrit( .1 * xmax, 1.1 * ymax )
          end if
c plot data
          do k = 1, kgrid
            do i = 1, nsteps
              d = 0
              do j = 1, ndimen
                d = d + sat( j + ( k - 1 ) * ndimen, i )**2
              end do
              d = sqrt( d )
              xs( i ) = ts * real( is( i ) )
              ys( i ) = d
            end do
            call jscoin( xs, ys, nsteps, krun )
          end do
        end if
      end if
c switch to next run or original if done
    4 krun = krun + 1
      if( krun .le. nruns )then
        call cmprun( krun )
      else
        call cmprun( 1 )
        krun = 1
      end if
      go to 5
    6 return
      end
