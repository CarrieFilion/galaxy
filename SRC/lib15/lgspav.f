      subroutine lgspav
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot time-averaged, rms & median values of Fourier coefficients
c   from a selected subset of the available data
c
c Data types handled by this routine are:
c      danl - sectorial harmonics of the mass distbn on grid rings
c      lgsp - logarithmic spiral coeffs
c      sfpc - coeffs saved from the SFP expansion
c      sphb - spherical Bessel functions coeffs
c      zanl - sectorial harmonics of the mid-plane displacements
c
c Called from ANALYS
c Graphics routines from JSPLOT
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
      include 'inc/bdpps.f'
c
      include 'inc/comprns.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      real grofu, roundup
c      real*8 gmassd, rofcx
c
c local allocatable arrays
      real, allocatable :: amax(:)
      real, allocatable :: amin(:)
      real, allocatable :: median(:)
      real, allocatable :: sigav(:)
      real, allocatable :: tav(:)
      real, allocatable :: tav2(:)
      real, allocatable :: sort(:)
c
c local array
      integer ifun( mostn + 1 )
c
c local variables
      character*4 type
      integer ifail, ip, it, j, k, krun, m, mact, mt, n
      logical first
      real a, prob, x, xmax, xmin, y, ymax, ymin
c      real f1
      real*8 r1, r2
c      real*8 r
      include 'inc/pi.f'
c
c select the disc population
      do ip = 1, ncmp
        if( disc( ip ) )icmp = ip
      end do
c select type of data
      call datype( type )
c
      call readin( type )
      print *, nt, ' moments in file - last is at time ', tme( nt )
c
c select Fourier harmonic
    2 call choosm( type, mact, m )
      if( mact .lt. 0 )return
c select data
      call select( mact, type )
c allocate space
      allocate ( amax( mr ) )
      allocate ( amin( mr ) )
      allocate ( median( mr ) )
      allocate ( sigav( mr ) )
      allocate ( tav( mr ) )
      allocate ( tav2( mr ) )
      allocate ( sort( nt ) )
      krun = 1
c find means, medians etc.
    5 if( krun .eq. 1 )then
        ymin = 1000
        ymax = 0
      end if
      do ip = jp, kp
        sigav( ip ) = 0
        tav( ip ) = 0
        tav2( ip ) = 0
        n = ip - kp
        do it = jt, kt
          n = n + kp - jp + 1
          if( mact .gt. 0 )then
            a = sqrt( sdata( 1, n )**2 + sdata( 2, n )**2 )
          else
            a = abs( sdata( 1, n ) )
          end if
          sort( it ) = a
          tav( ip ) = tav( ip ) + a
          tav2( ip ) = tav2( ip ) + a * a
          if( danl )sigav( ip ) = sigav( ip ) + sdata( 3, n )
        end do
        ifail = 0
        call srtmrg( sort( jt ), kt - jt + 1 )
        amax( ip ) = sort( kt )
        amin( ip ) = sort( jt )
        if( krun .eq. 1 )then
          ymax = max( ymax, amax( ip ) )
          ymin = min( ymin, amin( ip ) )
        end if
        median( ip ) = sort( ( jt + kt - 1 ) / 2 + 1 )
      end do
c set up frame
      if( krun .eq. 1 )then
        call jspage
c        call jschsz( .2 )
        call jssize( .1, .9, .1, .8 )
        if( sfpc )then
          xmax = real( ncp ) - .5
          xmin = -.5
        else if( danl )then
          xmax = grofu( real( nr( jgrid ) - 1 ) ) / lscale
          xmin = 0
        else if( lgsp )then
          xmax = 12.5
          xmin = -xmax
        else if( sphb )then
          xmax = ma * ( m + 1 )
          xmin = 1
        else if( zanl )then
          if( c3d )xmax = min( xm, ym ) / lscale
          if( p3a .or. p3d )xmax = rgrid( jgrid ) / lscale
          xmin = 0
        end if
c
        if( ymin .le. 0. )then
          ymin = -10.
        else
          ymin = log( ymin )
          ymin = -roundup( -ymin )
        end if
        ymax = max( log( ymax ), 0. )
        ymax = roundup( ymax )
        call jscale( xmin, xmax, ymin, ymax )
c draw axes
        if( sfpc )call jsaxis( 'x', 'Basis function', 1 )
        if( danl .or. zanl )call jsaxis( 'x', 'Radius', 1 )
        if( lgsp )call jsaxis( 'x', 'tan\\gamma', 1 )
        if( sphb )call jsaxis( 'x', 'Spherical Bessel function', 1 )
        call jsaxis( 'y', 'ln A', 1 )
c write header
        if( nruns .eq. 1 )then
          call jsbldt( 'Run no' )
          call jsbldi( irun, 5 )
        end if
        if( sphb )then
          call jsbldt( ' l  =' )
        else
          call jsbldt( ' m  =' )
        end if
        call jsbldi( mact, 3 )
        x = .8 * xmin + .2 * xmax
        y = 1.13 * ymax - .13 * ymin
        call jswrit( x, y )
        call jsbldt( 'Averaged over times' )
        call jsbldf( tme( jt ), 6, 1 )
        call jsbldt( ' to' )
        call jsbldf( tme( kt ), 6, 1 )
        x = .8 * xmin + .2 * xmax
        y = 1.05 * ymax - .05 * ymin
        call jswrit( xmin, y )
      end if
c evaluate averages and plot
      n = kt - jt + 1
      mt = 5
      if( nruns .gt. 1 )then
        mt = 1
        call pgsci( krun + 1 )
      end if
      do it = 1, mt
        if( it .eq. 2 )call jsdash( 0, 1, 0, 1 )
        if( it .eq. 4 )call jsdash( 2, 2, 2, 2 )
        if( it .eq. 5 )call jsdash( 0, 0, 0, 0 )
        do ip = jp, kp
          if( sfpc )x = ip - 1
          if( danl )x = grofu( real( ip - 1 ) ) / lscale
          if( lgsp )x = .25 * real( ip - ncp / 2 - 1 )
          if( zanl )x = ( real( ip ) - .5 ) / ( drfac * lscale )
          if( sphb )x = ip
          if( it .eq. 1 )then
            if( danl )sigav( ip ) = sigav( ip ) / real( n * na )
            tav( ip ) = tav( ip ) / real( n )
            y = tav( ip )
          end if
          if( it .eq. 2 )y = amin( ip )
          if( it .eq. 3 )y = amax( ip )
          if( it .eq. 4 )y = median( ip )
          if( it .eq. 5 )then
            tav2( ip ) = sqrt( tav2( ip ) / real( n ) )
            y = tav2( ip )
          end if
          y = log( max( y, 1.e-10 ) )
          if( ip .eq. jp )call jsmove( x, y )
          call jsline( x, y )
        end do
      end do
c draw expected amplitude
      if( krun .eq. 1 )then
        call jsdash( 3, 2, 0, 2 )
        if( danl )then
          if( .not. disc( 1 ) )call crash( 'LGSPAV',
     +                                           'Pop 1 is not a disk' )
          a = pi * real( nbod ) / ( fmass( 1 ) * cmpmas( 1 ) )
          x = max( real( jp ) - 1.5, 0. )
          r2 = grofu( x ) / lscale
          first = .true.
          do ip = jp, kp
            r1 = r2
            r2 = grofu( real( ip ) - .5 ) / lscale
            if( sigav( ip ) .gt. 0. )then
              prob = .25 * pi / ( a * ( r2**2 - r1**2 ) * sigav( ip ) )
              prob = .5 * log( prob )
              prob = min( prob, ymax )
              x = grofu( real( ip - 1 ) ) / lscale
              if( first )call jsmove( x, prob )
              call jsline( x, prob )
              first = .false.
            end if
          end do
        end if
        if( lgsp )then
          prob = pi / real( 4 * nbod )
          prob = .5 * log( prob )
          call pgsci( 1 )
          call jsmove( xmin, prob )
          call jsline( xmax, prob )
          call pgsci( krun + 1 )
        end if
      end if
      call jsdash( 0, 0, 0, 0 )
c label each curve
      if( nruns .gt. 1 )then
        x = .05 * xmax + .95 * xmin
        y = 1 - .06 * real( krun )
        y = y * ymax + ( 1. - y ) * ymin
        call jsmove( x, y )
        x = x + .05 * ( xmax - xmin )
        call jsline( x, y )
        call pgsci( 1 )
        call jsbldi( irun, 5 )
        x = x + .01 * ( xmax - xmin )
        y = y - .01 * ( ymax - ymin )
        call jswrit( x, y )
      end if
      call jsebuf
c switch to next run or original if done
      krun = krun + 1
      if( krun .le. nruns )then
        call cmprun( krun )
        call rewnd
        call nextrec( type, ifail )
        print *, 'starting run', krun, irun
        if( jt .gt. 1 )then
          do it = 2, jt
            call nextrec( type, ifail )
          end do
        end if
c read in selected data
        n = 0
        denf = 0
        do it = jt, kt
          call nextrec( type, ifail )
c take out com drift from m = 0 zanl data
c        if( zanl .and. ( m .eq. 0 ) )then
c          x = 0
c          z = 0
c          j = 1
c          do ip = 1, mr
c            k = j + mr * ( ma + 1 )
c            x = x + wres( k )
c            z = z + wres( j ) * wres( k )
c            j = j + ma + 1
c          end do
c          z = z / x
c          j = 1
c          do ip = 1, mr
c            wres( j ) = wres( j ) - z
c            j = j + ma + 1
c          end do
c        end if
c copy desired data to data array
          do ip = jp, kp
            if( danl .or. lgsp .or. vlgs )then
              j = ( ip - 1 ) * ma + m
              if( danl )j = j + 1
              if( vlgs .and. avlgs )j = j + 2 * ma * mr
              k = j + mr * ma
            else if( dnst )then
              j = ip
              k = j
            else if( s3dc )then
              j = ( ( m + 1 ) * ( m + 2 ) ) / 2
              j = 4 * ( s3ntm * ( ip - 1 ) + j - 1 ) + 1
              k = j + 1
            else if( sfpc )then
              if( ( it .gt. jt ) .and. ( wres( 1 ) .ne. maxr )
     +            .and. ( ip .eq. jp ) )print *,
     +                  'Warning - change of maxr from', maxr, ' to',
     +                            wres( 1 ), ' for time', it
              maxr = max( maxr, wres( 1 ) )
              j = 2 * ifun( ip )
              k = j + 1
            else if( scfc )then
              j = mact * ncp + 2 * ip - 1
              k = j + 1
            else if( sphb )then
              j = ip + ma * ( ( m * ( m + 1 ) ) / 2 )
              k = j + ma * ( ( mr + 1 ) * ( mr + 2 ) ) / 2
            else if( zanl )then
              j = ( ip - 1 ) * ( ma + 1 ) + m + 1
              k = j + mr * ( ma + 1 )
            end if
            n = n + 1
            if( vfld )then
c sum over angles - m=0 component of mean radial velocity only for now
              x = 0
              y = 0
              k = ( ip - 1 ) * nprop * ma - nprop
              do j = 1, ma
                k = k + nprop
                x = x + wres( k + 1 )
                y = y + wres( k + 1 ) * wres( k + 4 )
              end do
              if( x .gt. 0. )then
                sdata( 1, n ) = y / x
              else
                sdata( 1, n ) = 0
              end if
            else
              sdata( 1, n ) = wres( j )
              sdata( 2, n ) = wres( k )
            end if
            if( mact .eq. 0 )sdata( 2, n ) = 0
            denf = denf + sdata( 1, n )**2 + sdata( 2, n )**2
            if( danl )then
              j = j - m
              sdata( 3, n ) = wres( j )
            end if
          end do
        end do
        go to 5
      else
c finish up
        call cmprun( 1 )
        call pgsci( 1 )
        if( allocated( amax ) )deallocate ( amax )
        if( allocated( amin ) )deallocate ( amin )
        if( allocated( median ) )deallocate ( median )
        if( allocated( sort ) )deallocate ( sort )
        if( allocated( sigav ) )deallocate ( sigav )
        if( allocated( tav ) )deallocate ( tav )
        if( allocated( tav2 ) )deallocate ( tav2 )
        go to 2
      end if
      end
