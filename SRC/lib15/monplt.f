      subroutine monplt
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c part of the analysis software
c various options are provided to examine the time evolution of the
c   energies, angular momenta, etc of the particles that were selected
c   in the simulation to be monitored
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
      integer lnblnk
      real roundup
      real*8 Emin
c
c local arrays
      integer mmonit, mtime, jstep( 2 )
      parameter ( mmonit = 20000, mtime = 2000 )
      logical keep( mmonit )
      real ang( mmonit, 3 ), ene( mmonit, 3 )
      real mnsa( mtime, 2 ), mnse( mtime, 2 ), t( mtime )
      real x( mtime ), y( mtime )
c
c local variables
      character Et1*8, Et2*8, Lz1*12, Lz2*12, label*8
      integer i, ifail, iopt, iskip, itype, iw, j, jtime, k, n, nav
      integer nlive, ntime
      logical disk
      real amn, amx, berr, bint, dmn, dmx, emn, emx, f, rate, slerr
      real slope, t1, t2, xmax, xmin, xp, ymax, ymin, yp, y1, y2
      real*8 da2( 2 ), de2( 2 ), eav( 2, mtime )
c
      Et1 = 'E(000.0)'
      Et2 = 'E(000.0)'
      if( threed )then
        Lz1 = 'L(000.0)'
        Lz2 = 'L(000.0)'
      else
        Lz1 = 'L_{z}(000.0)'
        Lz2 = 'L_{z}(000.0)'
      end if
      icmp = 1
      disk = disc( icmp ) .or. rigidh
      call jssize( .25, .98, .1, .9 )
c choose plotting option
    3 call gtintg(
     +   'Enter plotting option (1-6), -1 to stop or 0 for help', iopt )
      if( iopt .lt. 0 )return
      if( iopt .eq. 0 )then
        print *, 'Options are:'
        print *, '1: comparison of values at two different times'
        print *, '2: differences of values at two different times'
        print *, '3: mn. sq. changes from start as functions of time'
        print *, '4: monitor a single particle - 2-D only'
        print *, '5: mean square changes from from moment to moment'
        print *, '6: mean energies over time of multiple mass species'
        go to 3
      end if
      if( iopt .gt. 6 )go to 3
      call rewnd
      call nextrec( 'MONI', ifail )
c initialize
    1 iskip = 1
      jstep( 1 ) = 0
      nlive = 0
      do i = 1, ncmp
        if( nsp( i ) .gt. 0 )nlive = nlive + 1
      end do
c
      if( iopt .lt. 3 )then
        call gtintgs(
     +         'Enter jstep & kstep or two zeros to stop', jstep, 2 )
        if( jstep( 1 ) + jstep( 2 ) .le. 0 )go to 3
        if( jstep( 2 ) .le. jstep( 1 ) )go to 1
      else if( iopt .eq. 4 )then
        call gtintg( 'Enter particle number to watch, 0 to end', iw )
        if( iw .le. 0 )go to 3
        if( iw .gt. ma )go to 1
        amn = 100
        amx = -amn
        dmn = 100
        dmx = -100
        emn = 100
        emx = -100
      else if( iopt .eq. 5 )then
        call gtintg( 'Enter stride through data, 1 for all', iskip )
      else if( iopt .eq. 6 )then
        if( nlive .le. 1 )then
          print *, 'not multiple populations', ncmp
          go to 3
        else if( nlive .gt. 2 )then
          call crash( 'MONPLT', 'Too many pops' )
        end if
      end if
c
      call rewnd
      ntime = 0
c read in data
    2 call nextrec( 'MONI', ifail )
      if( iskip .gt. 1 )then
        do i = 2, iskip
          if( ifail .eq. 0 )call nextrec( 'MONI', ifail )
        end do
      end if
c check space
      if( ma .gt. mmonit )call space( mmonit, ma, 'arrays', 'MONPLT' )
c mean square changes vs time
      if( ( iopt .eq. 3 ) .or. ( iopt .gt. 4 ) )then
        if( ifail .eq. 0 )then
          ntime = ntime + 1
          t( ntime ) = time
          if( ntime .eq. 1 )then
            do i = 1, ma
              j = mr * ( i - 1 )
              ene( i, 1 ) = wres( j + 1 )
              ang( i, 1 ) = wres( j + 2 )
              ang( i, 2 ) = wres( j + 3 )
              ang( i, 3 ) = wres( j + 4 )
              keep( i ) = ene( i, 1 ) .lt. 0.
            end do
            amx = 0
            emx = 0
            do k = 1, nlive
              mnsa( 1, k ) = 0
              mnse( 1, k ) = 0
            end do
            if( iopt .eq. 5 )then
              amn = 0
              emn = 0
              nav = 0
            end if
            if( iopt .ne. 6 )go to 2
          end if
c check space
          if( ntime .gt. mtime )call space(
     +                               mtime, ntime, 'arrays', 'MONPLT' )
c compute mean square changes
          n = 0
          do i = 1, 2
            da2( i ) = 0
            de2( i ) = 0
            eav( i, ntime ) = 0
          end do
          do i = 1, ma
            j = mr * ( i - 1 )
c            keep( i ) = keep( i ) .and.
c     +                  abs( wres( j + 1 ) - ene( i, 1 ) ) .lt. .5
c            if( keep( i ) )keep( i ) =
c     +   abs( ( wres( j + 2 ) - ang( i, 1 ) ) / ang( i, 1 ) ) .lt. .5
            if( c2d .and. keep( i ) )then
              f = sqrt( wres( j + 3 )**2 + wres( j + 4 )**2 )
              keep( i ) = f .lt. lscale
            end if
            if( ( wres( j + 1 ) .lt. 0. ) .and. keep( i ) )then
              n = n + 1
              k = 1
              if( nlive .gt. 1 )k = ( ( i - 1 ) * nlive ) / ma + 1
              if( ntime .gt. 1 )then
                de2( k ) = de2( k ) + ( wres( j + 1 ) - ene( i, 1 ) )**2
                da2( k ) = da2( k ) + ( wres( j + 2 ) - ang( i, 1 ) )**2
     +                              + ( wres( j + 3 ) - ang( i, 2 ) )**2
     +                              + ( wres( j + 4 ) - ang( i, 3 ) )**2
              end if
              eav( k, ntime ) = eav( k, ntime ) + wres( j + 1 )
            end if
            if( iopt .eq. 5 )then
              ene( i, 1 ) = wres( j + 1 )
              ang( i, 1 ) = wres( j + 2 )
              keep( i ) = wres( j + 1 ) .lt. 0.
            end if
          end do
          n = n / nlive
          do k = 1, nlive
            if( ntime .gt. 1 )then
              mnsa( ntime, k ) = sqrt( da2( k ) / dble( n ) )
c              mnsa( ntime, k ) = da2( k ) / dble( n )
              amx = max( amx, mnsa( ntime, k ) )
              mnse( ntime, k ) = de2( k ) / dble( n )
              emx = max( emx, mnse( ntime, k ) )
            end if
            eav( k, ntime ) = eav( k, ntime ) / dble( n )
            if( ( iopt .eq. 5 ) .and. ( time .gt. 2.01 ) )then
              amn = amn + mnsa( ntime, k )
              emn = emn + mnse( ntime, k )
              if( k .eq. 1 )nav = nav + 1
            end if
          end do
          go to 2
c end of file
        end if
        print *, ntime, ' values available'
        call gtintg( 'Enter initial number to omit from fit', jtime )
c work over both integrals
        do itype = 1, 2
c set scales and write title
          call jspage
          xmin = t( 1 )
          xmax = t( ntime )
          if( iopt .eq. 3 )then
            ymin = 0
            if( itype .eq. 1 )ymax = roundup( amx )
            if( itype .eq. 2 )ymax = roundup( emx )
            ymax = max( ymax, 1.e-5 )
c          if( idfn .eq. 11 )then
c            if( iopt .eq. 3 )ymax = 0.2
c            if( iopt .eq. 5 )ymax = 0.05
c          end if
          else
            ymin = 1000
            ymax = -ymin
            do i = 1, ntime
              do k = 1, nlive
                ymin = min( ymin, sngl( eav( k, i ) ) )
                ymax = max( ymax, sngl( eav( k, i ) ) )
              end do
            end do
            f = 2.1 * abs( ymax - ymin ) / abs( ymax + ymin )
            f = min( f, .1 )
            ymin = ymin * ( 1. + f )
            ymax = ymax * ( 1. - f )
          end if
          call jscale( xmin, xmax, ymin, ymax )
c write header
          if( itype .eq. 1 )then
            call jsbldt( 'Run no' )
            call jsbldi( irun, 4 )
            call jsbldt( '   N =' )
            f = nbod
            i = log10( f )
            j = nbod / 10**i
            label = ' x10^{ }'
            write( label( 1:1 ), '( i1 )' )j
            write( label( 7:7 ), '( i1 )' )i
            call jsbldt( label )
            if( nlive .gt. 1 )then
              call jsbldi( nlive, 2 )
              call jsbldt( 'species' )
            end if
            if( s3d )then
              call jsbldt( ' s3d grid l_{max} =' )
              call jsbldi( s3lmax, 2 )
            end if
            t1 = .9 * xmin + .1 * xmax
            t2 = 1.05 * ymax - .05 * ymin
            call jswrit( t1, t2 )
          end if
c plot results
          if( iopt .eq. 3 .or. iopt .eq. 5 )then
            call jsaxis( 'x', 'time', 1 )
            do k = 1, nlive
              if( itype .eq. 1 )then
                if( k .eq. 1 )call jsaxis( 'y', 'rms change in L', 1 )
                if( nlive .gt. 1 )call pgsci( k + 1 )
                call jsymbl( t, mnsa( 1, k ), ntime, 4 )
                call jsjoin( t, mnsa( 1, k ), ntime )
              else
                if( k .eq. 1 )call jsaxis( 'y', 'ms change in E', 1 )
                if( nlive .gt. 1 )call pgsci( k + 1 )
                call jsymbl( t, mnse( 1, k ), ntime, 4 )
                call jsjoin( t, mnse( 1, k ), ntime )
                if( nlive .gt. 1 )then
                  y1 = 1 - .07 * real( k )
                  y1 = ( 1. - y1 ) * ymin + y1 * ymax
                  call jsymbl( .1 * xmax + .9 * xmin, y1, 1, 4 )
                  call pgsci( 1 )
                  call jsbldt( '\\mu =' )
                  i = nint( 10. * cmpmas( k ) )
                  call jsbldi( i, 1 )
                  y1 = y1 - .01 * ymax
                  call jswrit( .15 * xmax + .85 * xmin, y1 )
                end if
              end if
              call pgsci( 1 )
              if( ( iopt .eq. 3 ) .and. ( ntime .gt. jtime + 1 ) )then
c                if( itype .eq. 1 )call linlsq( t( jtime + 1 ),
c     +   mnsa( jtime + 1, k ), ntime - jtime, slope, bint, slerr, berr )
                if( itype .eq. 1 )call linlsq( t,
     +   mnse( jtime + 1, k ), ntime - jtime, slope, bint, slerr, berr )
                if( itype .eq. 2 )call linlsq( t( jtime + 1 ),
     +   mnse( jtime + 1, k ), ntime - jtime, slope, bint, slerr, berr )
                rate = slope
                t1 = t( jtime + 1 )
                y1 = bint + slope * t1
                t2 = t( ntime )
                y2 = bint + slope * t2
                call jsdash( 2, 2, 2, 2 )
                call jsmove( t1, y1 )
                call jsline( t2, y2 )
                call jsdash( 0, 0, 0, 0 )
                call jsbldt( 'Rate from last' )
                call jsbldi( ntime - jtime, 3 )
                call jsbldt( 'values =' )
                call jsblde( rate, 8, 2 )
                t1 = .95 * xmin + .05 * xmax
                f = .05 * real( k )
                t2 = ( 1. - f ) * ymax + f * ymin
                call jswrit( t1, t2 )
              end if
            end do
            if( ( iopt .eq. 5 ) .and. ( nav .gt. 1 ) )then
              if( itype .eq. 1 )y1 = amn / real( nav )
              if( itype .eq. 2 )y1 = emn / real( nav )
              t1 = t( ntime + 1 - nav )
              t2 = t( ntime )
              call jsdash( 2, 2, 2, 2 )
              call jsmove( t1, y1 )
              call jsline( t2, y1 )
              call jsdash( 0, 0, 0, 0 )
              if( itype .eq. 1 )amn = amn / ( t2 - t1 )
              if( itype .eq. 2 )emn = emn / ( t2 - t1 )
              call jsbldt( 'Rate from last' )
              call jsbldi( nav, 3 )
              call jsbldt( 'values =' )
              if( itype .eq. 1 )call jsblde( amn, 8, 2 )
              if( itype .eq. 2 )call jsblde( emn, 8, 2 )
              t1 = .95 * xmin + .05 * xmax
              t2 = 1.05 * ymax - .05 * ymin
              call jswrit( t1, t2 )
            end if
          else if( iopt .eq. 6 )then
            call jsaxis( 'x', 'time', 1 )
            call jsaxis( 'y', '<E>', 1 )
            do k = 1, nlive
              call jsbldt( '\\mu =' )
              i = 10. * cmpmas( k )
              if( cmpmas( 3 ) .gt. 0. ) then
                i = nint( cmpmas( 1 ) / cmpmas( 2 ) )
                if( i .le. 9 )then
                  i = nint( cmpmas( k ) / cmpmas( 2 ) )
                  call jsbldi( i, 1 )
                else
                  i = nint( cmpmas( k ) / cmpmas( 2 ) )
                  call jsbldi( i, 2 )
                end if
              end if
              do i = 1, ntime
                y( i ) = eav( k, i )
              end do
              call linlsq( t, y, ntime, slope, bint, slerr, berr )
              if( nlive .gt. 1 )call pgsci( k + 1 )
              do i = 1, ntime
                f = eav( k, i )
                call jsymbl( t( i ), f, 1, 4 )
              end do
              if( nlive .gt. 1 )then
                f = 1. - .05 * real( k )
                yp = ( 1. - f ) * ymin + f * ymax
                xp = 4
                call jsymbl( xp, yp, 1, 4 )
                call pgsci( 1 )
c                call jsbldt( '\\mu =' )
                i = nint( 10. * cmpmas( k ) )
                call jsbldi( i, 1 )
                yp = yp - .01 * ( ymax - ymin )
                call jswrit( .15 * xmax + .85 * xmin, yp )
              end if
              call pgsci( 1 )
            end do
            go to 3
          end if
        end do
        go to 3
      end if
c single particle vs time
      if( iopt .eq. 4 )then
        if( ifail .eq. 0 )then
          ntime = ntime + 1
          t( ntime ) = time
c check space
          if( ntime .gt. mtime )call space(
     +                               mtime, ntime, 'arrays', 'MONPLT' )
c pick out selected particle
          j = mr * ( iw - 1 )
          mnsa( ntime, 1 ) = wres( j + 2 )
          amn = min( amn, mnsa( ntime, 1 ) )
          amx = max( amx, mnsa( ntime, 1 ) )
          mnse( ntime, 1 ) = wres( j + 1 )
          emn = min( emn, mnse( ntime, 1 ) )
          emx = max( emx, mnse( ntime, 1 ) )
          x( ntime ) = wres( j + 3 )
          y( ntime ) = wres( j + 4 )
          go to 2
        end if
c set scales and write title
        call jspage
        xmin = t( 1 )
        xmax = t( ntime )
        ymin = -roundup( -amn )
        ymax = roundup( amx )
        if( ymax .le. ymin )then
          ymin = roundup( amx - amn )
          ymax = .5 * ( amx + amn ) + ymin
          ymin = .5 * ( amx + amn ) - ymin
        end if
        call jscale( xmin, xmax, ymin, ymax )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        call jsbldt( 'particle no' )
        call jsbldi( iw, 4 )
        t1 = .9 * xmin + .1 * xmax
        t2 = 1.1 * ymax - .1 * ymin
        call jswrit( t1, t2 )
c plot results
        call jsaxis( 'x', 'time', 1 )
        call jsaxis( 'y', 'L_z', 1 )
        call jsymbl( t, mnsa, ntime, 4 )
        call jsjoin( t, mnsa, ntime )
c
        call jspage
        ymin = -roundup( -emn )
        ymax = roundup( emx )
        if( ymax .le. ymin )then
          ymin = roundup( emx - emn )
          ymax = .5 * ( emx + emn ) + ymin
          ymin = .5 * ( emx + emn ) - ymin
        end if
        call jscale( xmin, xmax, ymin, ymax )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        call jsbldt( 'particle no' )
        call jsbldi( iw, 4 )
        t1 = .9 * xmin + .1 * xmax
        t2 = 1.1 * ymax - .1 * ymin
        call jswrit( t1, t2 )
        call jsaxis( 'x', 'time', 1 )
        call jsaxis( 'y', 'E', 1 )
        call jsymbl( t, mnse, ntime, 4 )
        call jsjoin( t, mnse, ntime )
c
        call jspage
        t1 = rmax
        call jsescl( -t1, t1, -t1, t1 )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        call jsbldt( 'particle no' )
        call jsbldi( iw, 4 )
        call jswrit( -.9 * t1, 1.1 * t1 )
        call jsaxis( 'x', 'x', 1 )
        call jsaxis( 'y', 'y', 1 )
        call jsymbl( x, y, 1, 4 )
        call jsjoin( x, y, ntime )
c        call jspage
        go to 1
      end if
c store selected times only
      if( ( ifail .ne. 0 ) .or. ( istep .gt. jstep( 2 ) ) )then
        print *, 'Selected times not found'
        go to 1
      end if
c
      if( istep .eq. jstep( 1 ) )then
        t1 = time
        emn = 100
        emx = -100
        do i = 1, ma
          j = mr * ( i - 1 )
          if( wres( j + 1 ) .lt. 0. )then
            ene( i, 1 ) = wres( j + 1 )
            emn = min( emn, wres( j + 1 ) )
            emx = max( emx, wres( j + 1 ) )
            ang( i, 1 ) = wres( j + 2 )
            if( threed )then
              ang( i, 2 ) = wres( j + 3 )
              ang( i, 3 ) = wres( j + 4 )
            end if
          else
            ene( i, 1 ) = 0
            ang( i, 1 ) = 0
          end if
        end do
      end if
      if( istep .ne. jstep( 2 ) )go to 2
      t2 = time
      dmn = 0
      dmx = -1
      amn = 100
      amx = 0
      do i = 1, ma
        j = mr * ( i - 1 )
        if( wres( j + 1 ) .lt. 0. )then
          if( iopt .eq. 1 )then
            ene( i, 2 ) = wres( j + 1 )
            ang( i, 1 ) =
     +          sqrt( ang( i, 1 )**2 + ang( i, 2 )**2 + ang( i, 3 )**2 )
            if( threed )then
              ang( i, 2 ) =
     +    sqrt( wres( j + 2 )**2 + wres( j + 3 )**2 + wres( j + 4 )**2 )
            else
              ang( i, 2 ) = wres( j + 2 )
            end if
            amx = max( amx, ang( i, 2 ) )
            emn = min( emn, ene( i, 2 ) )
            emx = max( emx, ene( i, 2 ) )
          else
            ene( i, 3 ) = wres( j + 1 ) - ene( i, 1 )
            if( threed )then
              f = ( wres( j + 2 ) - ang( i, 1 ) )**2
     +          + ( wres( j + 3 ) - ang( i, 2 ) )**2
     +          + ( wres( j + 4 ) - ang( i, 3 ) )**2
              ang( i, 1 ) =
     +          sqrt( ang( i, 1 )**2 + ang( i, 2 )**2 + ang( i, 3 )**2 )
              ang( i, 3 ) = f
            else
              ang( i, 3 ) = wres( j + 2 ) - ang( i, 1 )
            end if
            amx = max( amx, ang( i, 3 ) )
            dmx = max( dmx, ene( i, 3 ) )
            amn = min( amx, ang( i, 3 ) )
            dmn = min( dmn, ene( i, 3 ) )
          end if
        else
          ene( i, 2 ) = 0
          ang( i, 2 ) = 0
          ene( i, 3 ) = 0
          ang( i, 3 ) = 0
        end if
      end do
c
      if( threed )then
        write( Lz1( 3:7 ), '( f5.1 )' )t1
        write( Lz2( 3:7 ), '( f5.1 )' )t2
      else
        write( Lz1( 7:11 ), '( f5.1 )' )t1
        write( Lz2( 7:11 ), '( f5.1 )' )t2
      end if
      write( Et1( 3:7 ), '( f5.1 )' )t1
      write( Et2( 3:7 ), '( f5.1 )' )t2
      call jspage
      if( iopt .eq. 1 )then
c set scaling
c        xmax = Lztmx
        xmax = amx
        if( disk )then
          xmin = Lzrmn( icmp )
        else
          xmin = 0
        end if
        ymin = xmin
        ymax = xmax
        call jsescl( xmin, xmax, ymin, ymax )
        call jsaxis( 'x', Lz1, 1 )
        call jsaxis( 'y', Lz2, 1 )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        t1 = .9 * xmin + .1 * xmax
        t2 = 1.1 * ymax - .1 * ymin
        call jswrit( t1, t2 )
c plot points
        n = ma / nlive
        j = 1
        do k = 1, nlive
          if( nlive .gt. 1 )call pgsci( k + 1 )
          call jsymbl( ang( j, 1 ), ang( j, 2 ), n, 0 )
          j = j + n
        end do
        call pgsci( 1 )
c
        call jspage
        xmin = Emin( 0.d0 )
c        xmax = emax( Lztmx )
        xmax = -1
        xmax = max( xmax, emx )
        xmin = emn
        ymin = xmin
        ymax = xmax
        call jsescl( xmin, xmax, ymin, ymax )
        call jsaxis( 'x', Et1, 1 )
        call jsaxis( 'y', Et2, 1 )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        t1 = .9 * xmin + .1 * xmax
        t2 = 1.1 * ymax - .1 * ymin
        call jswrit( t1, t2 )
c plot points
        n = ma / nlive
        j = 1
        do k = 1, nlive
          if( nlive .gt. 1 )call pgsci( k + 1 )
          call jsymbl( ene( j, 1 ), ene( j, 2 ), n, 0 )
          j = j + n
        end do
        call pgsci( 1 )
        go to 1
      end if
      if( iopt .eq. 2 )then
c set scaling
        if( disk )then
          xmin = Lzrmn( icmp )
        else
          xmin = 0
        end if
        xmax = Lztmx( icmp )
c        ymax = 0.1 * ( xmax - xmin )
        ymax = max( amx, -amn )
        ymin = -ymax
        if( threed )ymin = 0
        call jscale( xmin, xmax, ymin, ymax )
        call jsaxis( 'x', Lz1, 1 )
        k = lnblnk( Lz2 )
        call jsaxis( 'y', Lz2( 1:k )//'-'//Lz1, 1 )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        t1 = .9 * xmin + .1 * xmax
        t2 = 1.1 * ymax - .1 * ymin
        call jswrit( t1, t2 )
c plot points
        call jsymbl( ang, ang( 1, 3 ), ma, 0 )
c
        call jspage
c        xmin = Emin( 0.d0 )
c        xmax = Emax( Lztmx )
        xmax = -1
        xmax = max( xmax, emx )
        xmin = emn
c        ymax = 0.1 * ( xmax - xmin )
        ymax = max( dmx, -dmn )
        ymin = -ymax
        call jscale( xmin, xmax, ymin, ymax )
        call jsaxis( 'x', Et1, 1 )
        call jsaxis( 'y', Et2//'-'//Et1, 1 )
        call jsbldt( 'Run no' )
        call jsbldi( irun, 4 )
        t1 = .9 * xmin + .1 * xmax
        t2 = 1.1 * ymax - .1 * ymin
        call jswrit( t1, t2 )
c plot points
        call jsymbl( ene, ene( 1, 3 ), ma, 0 )
        go to 1
      end if
      go to 1
      end
