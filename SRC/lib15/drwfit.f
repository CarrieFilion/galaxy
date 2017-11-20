      subroutine drwfit( mact, nmodes, sumsq )
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c
c calling arguments
      integer mact, nmodes
      real sumsq
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
      include 'inc/grids.f'
c
c externals
      real grofu, roundup
c
c local variables
      integer imode, ip, iplot, it, n, nplot
      real ampm, cp, efac, p, phase, pm, sp, t, x, xmax
      real xmin, y, ymax, ymin
c
c set scaling
      ymax = 0.
      ymin = 100.
c find smallest data value
      do it = jt, kt
        n = ( it - jt ) * ( kp - jp + 1 )
        do ip = jp, kp
          n = n + 1
          y = sdata( 1, n )**2 + sdata( 2, n )**2
          if( y .gt. 0. )ymin = min( y, ymin )
        end do
      end do
      ymin = .5 * log( ymin )
      ymin = -roundup( -ymin )
c
      if( danl .or. dnst .or. s3dc .or. vfld )then
        xmax = rgrid( jgrid ) / lscale
        xmin = 0.
        if( mact .eq. 0 )ymax = -ymin
      else if( lgsp .or. vlgs )then
        xmax = 12.5
        xmin = -xmax
      else if( sfpc )then
        xmax = real( ncp ) - 0.5
        xmin = -0.5
        if( basset .eq. 'bess' )then
          xmax = real( ncp ) * deltak
          xmin = 0
        end if
        if( sfpc )ymin = -12
      else if( scfc .or. sphb )then
        xmax = real( ncp ) + 0.5
        xmin = 0.5
      else if( zanl )then
         xmin = .5 / ( drfac * lscale )
         xmax = ( real( ncp ) - .5 ) / ( drfac * lscale )
      end if
      if( xmin .ge. xmax )call crash( 'DRWFIT',
     +                                        'Unrecognized data type' )
c write heading
      call jspage
      call jssize( 0., .75, 0., 1. )
      call jscale( 0., 1., 0., 1. )
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      if( sphb )then
        call jsbldt( 'l =' )
      else
        call jsbldt( 'm =' )
      end if
      call jsbldi( mact, 2 )
      call jsbldt( 'Time interval' )
      x = ( tme( kt ) - tme( jt ) ) / real( kt - jt )
      call jsbldf( x, 6, 1 )
      call jswrit( .1, .93 )
      if( nmodes .eq. 0 )then
        nplot = 1
      else
        if( nmodes .eq. 1 )then
          nplot = 3
        else
          nplot = 4
        end if
        call jsbldi( nmodes, 2 )
        call jsbldt( 'Modes fitted to data' )
        if( nmodes .gt. 0 )then
          call jsbldt( 'Sum squares' )
          call jsbldf( sumsq, 7, 4 )
        end if
        call jswrit( .2, .87 )
      end if
c scan over plots
      do iplot = 1, nplot
c set up frame
        if( iplot .gt. 1 )call jspage
        call jssize( .1, .9, .1, .85 )
        call jscale( xmin, xmax, ymin, ymax )
c draw and mark axes
        if( danl .or. dnst .or. s3dc .or.
     +      vfld .or. zanl )call jsaxis( 'X', 'Radius', 1 )
        if( lgsp .or. vlgs )call jsaxis( 'x', 'tan\\gamma', 1 )
        if( sfpc .or. scfc )then
          if( basset .eq. 'bess' )then
            call jsaxis( 'x', 'k', 1 )
          else
            call jsaxis( 'x', 'Basis function', 1 )
          end if
        end if
        if( sphb )call jsaxis( 'x', 'Spherical Bessel function', 1 )
        call jsaxis( 'y', 'ln(A)', 1 )
c plot amplitude as a function of p
        do it = jt, kt
          n = ( it - jt ) * ( kp - jp + 1 )
          ampm = 0
          do ip = jp, kp
            if( danl .or. dnst.or. s3dc )
     +                          p = grofu( real( ip - 1 ) ) / lscale
            if( lgsp .or. vlgs )p = .25 * real( ip - ncp / 2 - 1 )
            if( vfld .or. zanl )
     +                   p = ( real( ip ) - .5 ) / ( drfac * lscale )
            if( sfpc )then
              p = ip - 1 + minn
              if( basset .eq. 'bess' )p = p * deltak
            end if
            if( sphb )p = ip + minn
            if( scfc )p = ip
            n = n + 1
c raw data
            if( iplot .eq. 1 )then
              x = sdata( 1, n )
              y = sdata( 2, n )
            end if
            if( ( iplot .eq. 2 ) .or. ( iplot .eq. 3 ) )then
c compute fitted value
              t = tme( it ) - tme( kt )
              x = 0
              y = 0
              do imode = 1, nmodes
                phase = freq( 1, imode ) * t
                cp = cos( phase )
                sp = sin( phase )
                efac = exp( freq( 2, imode ) * t )
                x = x + efac *
     +                 ( alp( ip, imode ) * cp - bet( ip, imode ) * sp )
                y = y + efac *
     +                 ( alp( ip, imode ) * sp + bet( ip, imode ) * cp )
              end do
              if( iplot .eq. 3 )then
                x = x - sdata( 1, n )
                y = y - sdata( 2, n )
              end if
            end if
c draw modes separately
            if( iplot .eq. 4 )then
              imode = it - jt + 1
              if( imode .gt. nmodes )return
              if( ip .eq. jp )ampm = 0.
              x = alp( ip, imode )
              y = bet( ip, imode )
              cp = sqrt( x * x + y * y )
              if( cp .gt. ampm )then
                ampm = cp
                pm = p
              end if
            end if
c plot value
            y = sqrt( x * x + y * y )
            y = max( y, .1e-10 )
            y = log( y )
            y = min( y, ymax )
c            y = max( y, ymin )
            if( ip .eq. jp )call jsmove( p, y )
            call jsline( p, y )
          end do
c mark maximum
          if( iplot .eq. 4 )then
            y = log( ampm )
            x = max( y, ymin )
            call jsymbl( pm, x, 1, 4 )
            call jsbldt( 'Mode no' )
            call jsbldi( imode, 3 )
c            call jsbldt( 'p = ' )
c            call jsbldf( pm, 4, 1 )
            call jswrit( p, y )
            if( y .ne. x )then
              call jsbldt( 'Max val' )
              call jsbldf( y, 6, 2 )
              call jswrit( p, y + .05 * ymin )
            end if
          end if
        end do
      end do
      return
      end
