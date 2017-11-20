      subroutine prpplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot velocity field data from the runXXX.res file
c    the data for disc particles have been binned differently from those
c    for halo particles
c There are 5 types of data for 2-D discs and 9 for 3-D discs or halos
c The surface density and Q profile plotting options work only for polar grids
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
      include 'inc/comprns.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      external gmassd, gmassh, lgsigma, sigmau, sigmav, qtoom, vcirc
      external vdisc, vhalo, vmean
      integer iscomp
      logical gtlogl
      real grofu, roundup, uofgr
      real*8 vcirc
c
c local arrays
      real, allocatable :: f( : )
      real, allocatable :: surf( : )
c
c local variables
      character type*4, xlabel*50, ylabel*50
      integer i, ia, iend, ifail, iopt, istrt, j, k, krun, l, ll, m, mh
      integer nc, nk, nrf, ns, ntime
      logical compr, compt, dtran, first, lknwn, OK
      real an, dr, dt, du, fnp, r, rmn, rsel, r1, r2, srf, u, vm, vm2
      real vo, tf, tm, tmax, tms, w, x, x1, x2, xmax, xmin, y, y1, y2
      real yf, ymax, ymin
      equivalence ( w, mh )
      include 'inc/pi.f'
c
c select active population
      call selpop( icmp )
      i = igrd( icmp )
      call switch( i )
      if( disc( icmp ) )then
        type = 'VFLD'
      else
        type = 'VFLH'
      end if
c determine whether theoretical curves are possible
      lknwn = .true.
      do i = 1, ncmp
        lknwn = lknwn .and. ( ctype( i ) .ne. 'UNKN' )
      end do
c
      if( nruns .gt. 1 )then
        compr = gtlogl( 'Compare time evolution at a fixed radius' )
        compt = .not. compr
c find duration of data
        ifail = 0
        call rewnd
        ia = 0
        tmax = 0
        do while ( ifail .eq. 0 )
          call nextrec( type, ifail )
          if( ifail .eq. 0 )ia = ia + 1
          tmax = max( tmax, time )
        end do
        dt = 0
        if( ia .gt. 1 )dt = tmax / real( ia - 1 )
        tmax = max( tmax, 1. )
        call rewnd
      else
        compt = .false.
        compr = .false.
      end if
c select best estimator of disc surface density
      call rewnd
      call nextrec( 'SIGR', ifail )
      sigr = ifail .eq. 0
      dtran = .not. sigr
      if( dtran )then
        call rewnd
        call nextrec( 'DANL', ifail )
        dtran = ifail .eq. 0
      end if
c
      ntime = 1
      if( nruns .eq. 1 )then
        call lstrid( type, ntime, iend )
      else
        iend = 10000000
        if( compr )then
          call gtreal( 'Enter radius to focus on', rsel )
          rsel = max( rsel, 1. / ( drfac * lscale ) )
        else
          OK = .false.
          do while ( .not. OK )
            call gtreal( 'Enter required analysis time', tms )
            OK = tms .eq. 0.
            i = 0
            if( dt .gt. 0. )then
              i = nint( tms / dt )
              an = tms - dt * real( i )
              OK = abs( an ) .lt. 0.01 * dt
            end if
            if( .not. OK )print *,
     +      'selected time is not available, nearest is', dt * real( i )
          end do
          tms = dt * real( i )
          print *, 'selected time', tms
        end if
      end if
c set scaling constants
      fnp = pmass / ( lscale**3 * ts**2 )
      rmn = rinh( jgrid )
c select option
    1 call jsebuf
      call gtintg( 'Enter option, -1 to stop or 0 for help', iopt )
      if( iopt .lt. 0 )return
c primitive help menu
      if( iopt .eq. 0 )then
        print *, 'Options are:'
        print *, ' 1 - m(r)'
        print *, ' 2 - <V> & V_c'
        print *, ' 3 - rms V dispersion'
        print *, ' 4 - <U>'
        print *, ' 5 - rms U dispersion'
        if( twod )then
          if( p2d )print *, ' 6 - Toomre Q'
          if( p2d )print *, ' 7 - surface density'
          print *, ' 8 - m(Lz)'
          print *, ' 9 - vortensity'
        end if
        if( threed )then
          print *, ' 6 - <z>'
          print *, ' 7 - rms z thickness'
          print *, ' 8 - <W>'
          print *, ' 9 - rms W dispersion'
          if( disc( icmp ) )then
            if( p3d )print *, '10 - Toomre Q'
            if( p3d )print *, '11 - surface density'
            print *, '12 - m(Lz)'
            print *, '13 - vortensity'
          else
            print *, '10 - m(Lz)'
          end if
        end if
        go to 1
      else if( ( iopt .gt. 13 ) .or. ( twod .and. iopt .gt. 9 ) )then
        go to 1
      else if( twod .and. ( iopt .gt. 5 ) )then
        iopt = iopt + 4
      else if( ( .not. disc( icmp ) ) .and. ( iopt .gt. 9 ) )then
        iopt = iopt + 2
      end if
c restart file
      krun = 1
    5 call rewnd
      nc = 0
c compute desired step for this run, allowing for different ts
      if( ( nruns .gt. 1 ) .and. compt )then
        ntime = nint( tms / ts )
        iend = ntime
      end if
      if( krun .eq. 1 )then
c initialise picture
        call jspage
c set window
        xmin = 0
        if( ( nruns .eq. 1 ) .or. compt )then
c          xmin = rmn / lscale
          xmax = rgrid( jgrid ) / lscale
c expand scale for condensed components
          xmax = min( xmax, 2. * rtrunc( icmp ) )
        else
          xmax = tmax
        end if
c
        ymax = 1.
        if( ( iopt .eq. 1 ) .or. ( iopt .eq. 12 ) )then
          ymax = real( nsp( icmp ) ) * fnp
          if( uqmass )ymax = ymax * fmass( icmp ) * cmpmas( icmp )
          ymax = max( ymax, 0.1 )
        else if( iopt .eq. 2 )then
          ymax = 2.
        else if( ( iopt .eq. 4 ) .or. ( iopt .eq. 8 ) )then
          ymax = .1
        else if( ( iopt .eq. 6 ) .or. ( iopt .eq. 7 ) )then
          ymax = roundup( .05 * rgrid( jgrid ) / lscale )
        else if( iopt .eq. 9 )then
          ymax = .5
        else if( iopt .eq. 13 )then
          ymax = .2
        end if
c double scaling for halo particles
        if( ( iopt .gt. 2  ) .and.
     +      ( iopt .lt. 10 ) .and.
     +      ( .not. disc( icmp ) ) )ymax = 2 * ymax
c scale up by "typical" velocity
        if( lknwn .and. ( iopt .gt. 1 ) .and. ( iopt .lt. 10 ) )then
          r = .5 * rtrunc( icmp )
          ymax = ymax * vcirc( dble( r ) )
        else if( iopt .eq. 10 )then
          ymax = 5.
        else if( iopt .eq. 11 )then
          ymax = 0.
        end if
c
        ymin = 0.
        if( ( iopt .eq. 4 ) .or. (iopt .eq. 6 ) .or.
     +      ( iopt .eq. 8 ) )then
          ymin = -ymax
        else if( iopt .eq. 11 )then
          ymin = -10.
        else if( iopt .eq. 12 )then
          call nextrec( 'ANGM', ifail )
          call rewnd
          xmax = real( mr - 2 ) / ( Lzfac * lscale**2 * ts )
          if( disc( icmp ) )then
            xmin = 0.
          else
            xmin = -xmax
          end if
        end if
        call jssize( .15, .95, .1, .9 )
        if( scale_set )then
          x1 = xmin
          x2 = xmax
          y1 = ymin
          y2 = ymax
          if( iopt .lt. 11 .or. iopt .eq. 13 )then
            x1 = x1 * unit_L
            x2 = x2 * unit_L
          end if
          if( iopt .eq. 1 )then
            y1 = y1 * unit_M
            y2 = y2 * unit_M
          else if( iopt .lt. 6 .or. iopt .eq. 8 .or. iopt .eq. 9 )then
            y1 = y1 * unit_V
            y2 = y2 * unit_V
          else if( iopt .eq. 6 .or. iopt .eq. 7 )then
            y1 = y1 * unit_L
            y2 = y2 * unit_L
          end if
          call jscale( x1, x2, y1, y2 )
        else
          call jscale( xmin, xmax, ymin, ymax )
        end if
c set y-axis label
        if( iopt .eq. 1 )then
          write( ylabel, '( a )' )'M(r)'
          i = 4
        else if( iopt .eq. 2  )then
          write( ylabel, '( a )' )'<V> & V_{c}'
          i = 11
        else if( iopt .eq. 3  )then
          write( ylabel, '( a )' )'\\sigma_{v}'
          i = 11
        else if( iopt .eq. 4  )then
          write( ylabel, '( a )' )'<U>'
          i = 3
        else if( iopt .eq. 5  )then
          write( ylabel, '( a )' )'\\sigma_{u}'
          i = 11
        else if( iopt .eq. 6  )then
          write( ylabel, '( a )' )'<z>'
          i = 3
        else if( iopt .eq. 7  )then
          write( ylabel, '( a )' )'\\sigma_{z}'
          i = 11
        else if( iopt .eq. 8  )then
          write( ylabel, '( a )' )'<W>'
          i = 3
        else if( iopt .eq. 9  )then
          write( ylabel, '( a )' )'\\sigma_{w}'
          i = 11
        else if( iopt .eq. 10 )then
          write( ylabel, '( a )' )'Q'
          i = 1
        else if( iopt .eq. 11 )then
          write( ylabel, '( a )' )'ln(\\Sigma)'
          i = 11
        else if( iopt .eq. 12 )then
          write( ylabel, '( a )' )'L_{z}'
          i = 5
        else if( iopt .eq. 13 )then
          write( ylabel, '( a )' )'vortensity'
          i = 10
        end if
c set x-axis label
        if( scale_set )then
          if( compt )then
            write( xlabel, '( a )' )'time (Myr)'
            j = 10
          else
            if( iopt .eq. 12 )then
              write( xlabel, '( a )' )'L_{z}'
              j = 5
            else
              write( xlabel, '( a )' )'R (kpc)'
              j = 7
            end if
          end if
        else
          if( compt )then
            write( xlabel, '( a )' )'time'
            j = 4
          else
            if( iopt .eq. 12 )then
              write( xlabel, '( a )' )'L_{z}'
              j = 5
            else
              write( xlabel, '( a )' )'R'
              j = 1
            end if
          end if
        end if
c draw axes
        call jsthik( 3 )
        if( compr )then
          if( iopt .eq. 1 )i = i - 3
          write( ylabel( i+1:i+5 ), '( a )' )'(   )'
          write( ylabel( i+2:i+4 ), '( f3.1 )' )rsel
          i = i + 5
        end if
        call jsaxis( 'X', xlabel( 1:j ), 1 )
        call jsaxis( 'Y', ylabel( 1:i ), 1 )
        call jscale( xmin, xmax, ymin, ymax )
      end if
c
c start of main loop
c
      i = krun + 1
      if( krun .gt. 8 )i = krun - 7
      if( nruns .gt. 1 )call pgsci( i )
      if( i .lt. krun )call jsdash( 2, 2, 2, 2 )
   20 continue
      if( nc .eq. 1 .and. krun .eq. 1 )istrt = istep
c density distribution
      if( ( iopt .eq. 10 ) .or. ( iopt .eq. 11 ) .or.
     +    ( iopt .eq. 13 ) )then
        if( sigr )then
          call nextrec( 'SIGR', ifail )
        else if( dtran )then
          call nextrec( 'DANL', ifail )
        else
          if( iopt .eq. 10 )go to 30
          call nextrec( 'DNST', ifail )
        end if
        if( ( ifail .ne. 0 ) .or. ( istep .gt. iend ) )go to 70
        if( compt )then
          if( istep .lt. ntime )go to 20
        else
          if( mod( istep, ntime ) .ne. 0 )go to 20
        end if
        if( iopt .eq. 11 )nc = nc + 1
        m = 0
        if( dtran )m = 1 - ma
c find surface density at selected radius
        if( sigr .and. compr .and. ( iopt .eq. 11 ) )then
          x = -rsel
          y = -1
          m = -1
c find straddling data
          do while( x * y .gt. 0. )
            m = m + 2
            y = x
            if( m .lt. mr )then
              x = wres( m ) - rsel
              if( m .gt. 1 )dr = wres( m ) - wres( m - 2 )
            else
              x = 1
            end if
          end do
c interpolate
          if( m .eq. 1 )then
            an = wres( 2 )
          else if( m .lt. mr )then
            x = x / dr
            an = x * wres( m - 1 ) + ( 1. - x ) * wres( m + 1 )
          else
            an = 0
          end if
c plot
          if( an .le. 0. )then
            y = ymin
          else
            y = max( ymin, log( an ) )
          end if
          if( istep .eq. 0 )call jsmove( time, y )
          call jsline( time, y )
        else
c work over all radii
          first = .true.
          r = 0
          ns = mr
          if( sigr )ns = nr( igrd( icmp ) )
          if( .not. allocated( surf ) )allocate ( surf( ns ) )
          do i = 1, ns
            an = 0.
            dr = r
            r = grofu( real( i - 1 ) ) / lscale
            if( sigr )then
c find straddling radii
              if( i .eq. 1 )r = 1.e-6
              x = -r
              y = -1
              m = -1
              do while( x * y .gt. 0. )
                m = m + 2
                y = x
                if( m .lt. mr )then
                  x = wres( m ) - r
                  if( m .gt. 1 )dr = wres( m ) - wres( m - 2 )
                else
                  x = 1
                end if
              end do
c interpolate
              if( m .eq. 1 )then
                an = wres( 2 )
              else if( m .lt. mr )then
                x = x / dr
                an = x * wres( m - 1 ) + ( 1. - x ) * wres( m + 1 )
              else
                an = 0
              end if
            else if( dtran )then
c extract monopole term
              dr = r - dr
              m = m + ma
              an = wres( m ) / real( na )
            else
c average around a ring
              if( p2d .or. p3d )then
                do ia = 1, na
                  m = m + 1
                  an = an + wres( m )
                end do
                an = an / real( na )
              end if
c no need to average
              if( p3a )then
                m = m + 1
                an = wres( m )
              end if
            end if
c save value
            surf( i ) = an
c plot
            if( iopt .eq. 11 )then
              if( an .le. 0. )then
                y = ymin
              else
                y = max( ymin, log( an ) )
              end if
              if( ( nruns .eq. 1 ) .or. compt )then
                if( first )call jsmove( r, y )
                first = .false.
                if( y .gt. ymin )then
                  call jsline( r, y )
                else
                  call jsmove( r, y )
                end if
              else if( abs( r - rsel ) .lt. dr )then
                if( istep .eq. 0 )call jsmove( time, y )
                call jsline( time, y )
              end if
            end if
          end do
        end if
        if( ( iopt .eq. 10 ) .or. ( iopt .eq. 13 ) )go to 30
        if( nc .eq. 1 .and. nruns .eq. 1 .and. lknwn )then
c draw theoretical curve
          call jsdash( 2, 2, 2, 2 )
          call jsplt2( lgsigma )
          call jsdash( 0, 0, 0, 0 )
        end if
        go to 20
      end if
c frequencies
      if( iopt .eq. 2 )then
   30   call nextrec( 'FRQS', ifail )
        if( ( ifail .ne. 0 ) .or. ( istep .gt. iend ) )go to 70
        if( compt )then
          if( istep .ne. ntime )go to 30
        else
          if( mod( istep, ntime ) .ne. 0 )go to 30
        end if
        if( .not. allocated( f ) )allocate( f( mr ) )
        first = .true.
        r = 0
        nrf = min( mr, size( f ) )
        do i = 1, nrf
c circular velocity
          if( iopt .eq. 2 )then
            dr = r
            r = ( rmn + real( i ) * drfrqs ) / lscale
            dr = r - dr
            y = r * wres( i )
            if( ( nruns .eq. 1 ) .or. compt )then
              if( first )call jsmove( r, y )
              first = .false.
              if( y .gt. 0. )then
                call jsline( r, y )
              else
                call jsmove( r, y )
              end if
            else if( abs( r - rsel ) .lt. .5 * dr )then
              if( first )then
                tf = time
                yf = y
                first = .false.
              end if
              call jsmove( tf, yf )
              call jsline( time, y )
              tf = time
              yf = y
            end if
          end if
c remember epicycle frequency for plotting Q later
          if( iopt .eq. 10 )then
            k = i + mr
            f( i ) = wres( k )
          end if
c vortensity
          if( iopt .eq. 13 )then
            r = ( rmn + real( i ) * drfrqs ) / lscale
            if( sigr .or. dtran )then
              u = uofgr( r * lscale )
              j = u
              du = u - real( j )
              srf = surf( j + 1 ) * ( 1. - du ) + surf( j + 2 ) * du
            else
              r1 = rmn + real( i - 1 ) / ( lscale * drfac )
              r2 = rmn + real( i ) / ( lscale * drfac )
              srf = an * fnp / ( pi * ( r2**2 - r1**2 ) )
            end if
            if( wres( i + mr ) .ne. 0. )then
              y = 2 * wres( i ) * srf / wres( i + mr )**2
            else
              y = 0
              first = .true.
            end if
            if( first )call jsmove( r, y )
            first = .false.
            if( y .gt. 0. )then
              call jsline( r, y )
            else
              call jsmove( r, y )
            end if
          end if
        end do
        nk = mr - 1
      end if
c angular momentum
      if( iopt .eq. 12 )then
        call nextrec( 'ANGM', ifail )
        if( ( ifail .ne. 0 ) .or. ( istep .gt. iend ) )go to 70
        if( compt )then
          if( istep .ne. ntime )go to 20
        else
          if( mod( istep, ntime ) .ne. 0 )go to 20
        end if
        nc = nc + 1
c set bin width & offset
        an = xmax / real( mr - 1 )
        if( .not. disc( icmp ) )an = 2 * an
c plot distribution
        x = xmin - an
        y = 0
        call jsmove( x, y )
        do i = 1, mr
          w = wres( i )
          y = y + fnp * real( mh )
          call jsline( x, y )
          x = x + an
          call jsline( x, y )
        end do
        go to 20
      end if
c velocity field data
   50 call nextrec( type, ifail )
      if( ( ifail .ne. 0 ) .or. ( istep .gt. iend ) )go to 70
      if( compt )then
        if( istep .ne. ntime )go to 50
      else
        if( mod( istep, ntime ) .ne. 0 )go to 50
      end if
      nc = nc + 1
c work over radii
      y = 0.
      first = .true.
      ia = 1
      if( iopt .eq. 10 )ia = 2
      r = 0
      do i = ia, mr
        dr = r
        r = rmn + ( real( i ) - .5 ) / drfac
        if( iopt .eq. 1 )r = rmn + real( i ) / drfac
        r = r / lscale
        dr = r - dr
        if( r .lt. xmax )then
c sum over angles or planes
          vm = 0.
          vm2 = 0.
          an = 0.
          l = ( i - 1 ) * nprop * ma - nprop
          do j = 1, ma
            l = l + nprop
            an = an + wres( l + 1 )
            ll = l + iopt
            if( iopt .eq. 10 )ll = l + 5
            if( ( iopt .eq. 2 ) .or. ( iopt .eq. 4 ) .or.
     +          ( iopt .eq. 6 ) .or. ( iopt .eq. 8 ) )then
              vm = vm + wres( ll ) * wres( l + 1 )
            else if( ( iopt .eq. 3 ) .or. ( iopt .eq. 5 ) .or.
     +               ( iopt .eq. 7 ) .or. ( iopt .eq. 9 ) .or.
     +               ( iopt .eq. 10 ) )then
              y = wres( ll - 1 )
              vm = vm + y * wres( l + 1 )
              vm2 = vm2 +
     +               wres( l + 1 ) * ( wres( ll ) * wres( ll ) + y * y )
            end if
          end do
          if( iopt .eq. 1 )then
            y = y + an * fnp
          else
c compute mean and dispersion for annulus
            y = 0.
            if( an .gt. 0. )then
              y = vm / an
              if( iopt .eq. 4 )then
                y = min( y, 1.2 * ymax )
                y = max( y, 1.2 * ymin )
              end if
              if( ( mod( iopt, 2 ) .ne. 0 ) .or. ( iopt .eq. 10 ) )then
                y = ( vm2 - an * y * y ) / an
                y = max( 0., y )
                y = sqrt( y )
              end if
c compute Q
              if( iopt .eq. 10 )then
                k = ( r * lscale - rmn ) / drfrqs
                k = max( k + 1, 2 )
                k = min( k, nk + 1 )
                if( sigr .or. dtran )then
                  u = uofgr( r * lscale )
                  j = u
                  if( j .lt. ns - 1 )then
                    du = u - real( j )
                    srf = surf( j + 1 ) * ( 1. - du ) +
     +                    surf( j + 2 ) * du
                  else
                    srf = 0
                  end if
                else
                  r1 = rmn + real( i - 1 ) / ( lscale * drfac )
                  r2 = rmn + real( i ) / ( lscale * drfac )
                  srf = an * fnp / ( pi * ( r2**2 - r1**2 ) )
                end if
                if( srf .eq. 0. )then
                  y = 0.
                else
c expression for Q
                  y = y * f( k ) / ( 3.36 * srf )
c expression for Lcrit
c                  y = 4 * pi**2 * srf / f( k )**2
                end if
                y = min( y, ymax )
c expression for mean epicylce size
c              else if( iopt .eq. 5 )then
c                k = ( r * lscale - rmn ) / drfrqs
c                k = max( k + 1, 2 )
c                k = min( k, nk + 1 )
c                y = y / f( k )
              end if
            end if
          end if
          if( ( nruns .eq. 1 ) .or. compt )then
c plot radial variation at selected times
            if( first )call jsmove( r, y )
            if( an .gt. 0. )then
              first = .false.
              call jsline( r, y )
            else
              call jsmove( r, y )
            end if
c plot time evolution at selected radius
          else if( abs( r - rsel ) .lt. .6 * dr )then
            if( nc .eq. 1 )then
              tm = time
              vo = y
            end if
            if( time .gt. tm )then
              call jsmove( tm, vo )
              call jsline( time, y )
              tm = time
              vo = y
            end if
          end if
        end if
      end do
c draw theoretical curves
      if( ( nc .eq. 1 ) .and. .not. ( ( iopt .eq. 4 ) .or.
     +         ( iopt .eq. 6 ) .or. ( iopt .eq. 8 ) ) .and.
     +         ( nruns .eq. 1 ) .and. lknwn )then
        call jsdash( 2, 2, 2, 2 )
        call pgsci( 1 )
        if( iopt .eq. 1 )then
          if( disc( icmp ) )then
            call jsplt2( gmassd )
          else
            call jsplt2( gmassh )
          end if
        end if
        if( iopt .eq. 2 )then
          call jsplt2( vcirc )
          if( ( ncmp .gt. 1 ) )then
            j = icmp
            k = icomp
            do icmp = 1, ncmp
              if( cmpmas( icmp ) .gt. 0. )then
                if( disc( icmp ) )then
                  call jsdash( 0, 1, 0, 1 )
                  call jsplt2( vdisc )
                else
                  if( cdft( icmp ) .eq. 'COMP' )icomp = iscomp( icmp )
                  call jsdash( 3, 1, 0, 1 )
                  call jsplt2( vhalo )
                end if
              end if
            end do
            icmp = j
            icomp = k
          end if
        end if
        if( dist( icmp ) .and. ( ( iopt .eq. 2 ) .or.
     +                ( iopt .eq. 5 ) .or. ( iopt .eq. 10 ) ) )then
          if( gtlogl( 'Plot theoretical curve?' ) )then
            if( iopt .eq. 2 )call jsplt2( vmean )
            if( iopt .eq. 5 )call jsplt2( sigmau )
            if( iopt .eq. 10 )call jsplt2( qtoom )
          end if
        end if
        call jsdash( 0, 0, 0, 0 )
      end if
      go to 20
c end of file
   70 call jsdash( 0, 0, 0, 0 )
      if( nruns .eq. 1 )then
        call jsbldt( 'Run no' )
        call jsbldi( irun, 5 )
        call jsbldt( 'time interval' )
        i = ntime
        if( ( ntime .eq. 1 ) .and.
     +      ( nc .gt. 1 ) )i = ( istep - istrt ) / ( nc - 1 )
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
      else
        x = .05 * xmax + .95 * xmin
        y = .04 * real( krun )
        y = y * ymax + ( 1. - y ) * ymin
        call jsmove( x, y )
        x = x + .05 * ( xmax - xmin )
        call jsline( x, y )
        call pgsci( 1 )
        call jsbldi( irun, 5 )
        x = x + .01 * ( xmax - xmin )
        y = y - .01 * ( ymax - ymin )
        if( ( krun .eq. 1 ) .and. compt  )then
          call jsbldt( 't=' )
          tm = ts * real( ntime )
          i = nint( tm )
          if( i .gt. 9999 )then
            call jsbldi( i, 5 )
          else if( i .gt. 999 )then
            call jsbldi( i, 4 )
          else if( i .gt. 99 )then
            call jsbldi( i, 3 )
          else
            call jsbldi( i, 2 )
          end if
        end if
        call jswrit( x, y )
      end if
c switch to the next run or the original if done
      krun = krun + 1
      if( krun .le. nruns )then
        call cmprun( krun )
        go to 5
      else
        krun = 1
        call cmprun( 1 )
        go to 1
      end if
      end
