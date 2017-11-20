      subroutine loganl
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot results from several different types of analysis saved
c   in the runXXX.res file.  The variation at specific times is plotted
c
c Data types handled by this routine are:
c      danl - sectorial harmonics of the mass distbn on grid rings
c      dnst - radial mass distribution on grid rings
c      lgsp - logarithmic spiral coeffs
c      s3dc - coeffs saved from the S3D grid
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
      real grofu, roundup, uofgr
c
c local arrays
      integer ifun( mostn + 1 )
      integer jr, mstp
      parameter ( jr = 201, mstp = 5001 )
      real amp( mstp ), phase( mstp ), tm( mstp )
      real ph( jr, mruns )
c
c local variables
      character type*4
      integer i, iend, ifail, ip, istp, j, jend, jstart, jstep, k, krun
      integer m, mact, man, np
      logical compt, allowed
      real a, p, x, xmax, xmin, y, ymax, ymin
      include 'inc/pi.f'
c
      data jend, jstart / 2 * 0 /
c
c select type of data
      call datype( type )
      if( type( 1:3 ) .eq. 'END' )return
c select Fourier harmonic
    1 call choosm( type, mact, m )
      if( mact .lt. 0 )return
      if( danl )m = m + 1
      if( sfpc )then
        np = 0
        do ip = 1, lastf
          if( mact .eq. msel( ip ) )then
            np = np + 1
            if( np .gt. mostn + 1 )call crash( 'LOGANL',
     +                                          'ifun array too small' )
            ifun( np ) = ip
          end if
        end do
        if( np .eq. 0 )call crash( 'LOGANL',
     +                        'No basis function with this m selected' )
      end if
c
      compt = ( nruns .gt. 1 )
    2 if( compt )then
        call gtreal( 'Enter required analysis time, or -ve to end', a )
        if( a .lt. 0 )go to 1
        iend = a / ts
        jstep = iend
      else
        call gtintg(
     +             'Enter interval between plots in time-steps', jstep )
        jstep = max( 1, jstep )
        call gtintg(
     +     'Enter step no of last moment required, or 0 for all', iend )
        if( iend .eq. 0 )iend = 1000000
      end if
      istp = 0
c set scaling
      ymax = 0.
      if( quiet( 1 ) .or. quiet( 2 ) )then
        ymin = -15.
      else
        ymin = -6.
      end if
      if( danl .or. dnst .or. s3dc )then
        xmin = 0
        xmax = rgrid( jgrid ) / lscale
      else if( lgsp )then
        xmax = 12.5
        xmin = -xmax
      else if( scfc .or. sfpc )then
        if( scfc )np = ncp
        xmax = real( np ) - .5
        xmin = -.5
        ymax = 0
      else if( sphb )then
        xmin = 0
        ymin = -10
        xmax = ma * ( m + 1 ) + 1
c work over radial orders and m at this l
        do ip = 1, ma * ( m + 1 )
          j = ip + ma * ( ( m * ( m + 1 ) ) / 2 )
          k = j + ma * ( ( mr + 1 ) * ( mr + 2 ) ) / 2
          a = sqrt( wres( j )**2 + wres( k )**2 )
          if( a .gt. 0. )ymin = min( ymin, log( a ) )
        end do
        ymin = -roundup( -ymin )
      else if( zanl )then
        xmin = 0.
        xmax = rgrid( jgrid ) / lscale
        ymin = -10
        do ip = 1, mr
          j = ( ip - 1 ) * ( ma + 1 ) + m + 1
          k = j + mr * ( ma + 1 )
          if( m .gt. 0 )then
            a = sqrt( wres( j ) * wres( j ) + wres( k ) * wres( k ) )
          else
            a = abs( wres( j ) )
          end if
          if( a .gt. 0. )ymin = min( ymin, log( a ) )
        end do
        if( ymin .lt. 0. )ymin = -roundup( -ymin )
      end if
      if( xmin .ge. xmax )call crash( 'LOGANL',
     +                                        'Unrecognized data type' )
      ymin = max( ymin, -20. )
c draw and mark axes
      call jspage
      call jssize( .15, .7, .1, .85 )
      call jscale( xmin, xmax, ymin, ymax )
      if( danl .or. dnst .or. zanl )call jsaxis( 'x', 'R', 1 )
      if( s3dc )call jsaxis( 'x', 'r', 1 )
      if( lgsp )call jsaxis( 'x', 'tan\\gamma', 1 )
      if( sfpc .or. scfc )call jsaxis( 'x', 'Basis function', 1 )
      if( sphb )call jsaxis( 'x', 'radial function', 1 )
      call jsaxis( 'y', 'ln(A)', 1 )
c work through file
      krun = 1
    5 if( nruns .gt. 1 )call pgsci( krun + 1 )
      ifail = 0
      do while ( ifail .eq. 0 )
        istp = istp + 1
        if( ( compt .and. ( istep .eq. iend ) ) .or.
     +      ( nruns .eq. 1 ) )then
c plot amplitude variation at this moment
          call jsthik( 1 )
          if( .not. ( scfc .or. sfpc ) )np = mr
          if( sphb )np = ma * ( m + 1 )
          do ip = 1, np
            if( dnst .or. danl .or. lgsp  )then
              if( danl .or. dnst )p = grofu( real( ip - 1 ) ) / lscale
              if( lgsp )p = .25 * real( ip - mr / 2 - 1 )
              j = ( ip - 1 ) * ma + m
              k = j + mr * ma
            else if( s3dc )then
              p = grofu( real( ip - 1 ) ) / lscale
              j = ( ( m + 1 ) * ( m + 2 ) ) / 2
              j = 4 * ( s3ntm * ( ip - 1 ) + j - 1 ) + 1
              k = j + 1
            else if( scfc )then
              p = ip - 1
              j = mact * np + 2 * ip - 1
              k = j + 1
            else if( sfpc )then
              p = ip - 1
              j = 2 * ifun( ip )
              k = j + 1
            else if( sphb )then
              p = ip
              j = ip + ma * ( ( m * ( m + 1 ) ) / 2 )
              k = j + ma * ( ( mr + 1 ) * ( mr + 2 ) ) / 2
            else if( zanl )then
              p = ( real( ip ) - .5 ) / ( drfac * lscale )
c              if( p .gt. xmax )go to 3
              j = ( ip - 1 ) * ( ma + 1 ) + m + 1
              k = j + mr * ( ma + 1 )
            end if
            if( m .gt. 0 )then
              a = sqrt( wres( j ) * wres( j ) + wres( k ) * wres( k ) )
              if( compt )then
                if( ip .gt. jr )call crash( 'LOGANL',
     +                                            'ph array too small' )
                ph( ip, krun ) =
     +                      atan2( wres( k ), wres( j ) ) / real( mact )
              end if
            else
              a = abs( wres( j ) )
            end if
c
            y = ymin
            if( a .gt. 0. )y = log( a )
            if( ip .eq. 1 )call jsmove( p, y )
            if( p .lt. xmax )call jsline( p, y )
          end do
        end if
c get next record needed record
        call nextrec( type, ifail )
        do while ( ( ifail .eq. 0 ) .and.
     +             ( mod( istep, jstep ) .ne. 0 ) )
          call nextrec( type, ifail )
        end do
        if( istep .gt. iend )ifail = 1
      end do
c write heading
      if( nruns .eq. 1 )then
        call jsbldt( 'Run no' )
        call jsbldi( irun, 5 )
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
        call jschsz( .25 )
        call jswrit( x, y )
        call jschsz( .3 )
      end if
      if( krun .eq. 1 )then
        if( sphb .or. s3dc .or. scfc )then
          call jsbldt( ' l = ' )
        else
          call jsbldt( ' m = ' )
        end if
        call jsbldi( mact, 2 )
        if( nruns .eq. 1 )then
          call jsbldt( ' time interval' )
          if( istp .gt. 1 )x = ts * real( istep ) / real( istp - 1 )
          call jsbldf( x, 6, 1 )
        else
          call jsbldt( 't=' )
          a = ts * real( iend )
          i = nint( a )
          call jsbldi( i, 3 )
        end if
        x = .8 * xmin + .2 * xmax
        y = -.05 * ymin + 1.05 * ymax
        call jsthik( 3 )
        call jswrit( x, y )
      end if
c plot phases at comparison time
      if( compt .and. ( krun .eq. nruns ) )then
        call jspage
        ymax = pi / real( mact )
        ymin = -ymax
        call jscale( xmin, xmax, ymin, ymax )
        call jsaxis( 'x', 'R', 1 )
        call jsaxis( 'y', '\\phi', 1 )
        do i = 1, nruns
          call pgsci( i + 1 )
          do ip = 1, np
            if( dnst .or. danl .or. lgsp  )then
              if( danl .or. dnst )p = grofu( real( ip - 1 ) ) / lscale
              if( lgsp )p = .25 * real( ip - mr / 2 - 1 )
            else if( s3dc )then
              p = grofu( real( ip - 1 ) ) / lscale
            else if( scfc .or. sfpc )then
              p = ip - 1
            else if( sphb )then
              p = ip
            else if( zanl )then
              p = ( real( ip ) - .5 ) / ( drfac * lscale )
            end if
            y = ph( ip, i )
            if( ip .eq. 1 )call jsmove( p, y )
            call jsline( p, y )
          end do
        end do
        call pgsci( 1 )
      end if
c select p value for growth rate
      allowed = ( mact .ne. 0 ) .and. ( nruns .eq. 1 )
      do while ( allowed )
        if( danl .or. s3dc .or. zanl )call gtreal( 'Select radius', p )
        if( lgsp )call gtreal(
     +       'Select tan(gamma) value or large value to skip', p )
        if( scfc .or. sfpc .or. sphb )then
          call gtreal( 'Select function', p )
          p = p + 1
        end if
        if( lgsp )then
          allowed = int( abs( p ) ) .lt. np
        else
          allowed = p .gt. 0.
        end if
        if( allowed )then
c rewind file and set pointers to pick out requested values
          call rewnd
          istp = 0
          call nextrec( type, ifail )
          if( danl )then
            ip = uofgr( p * lscale ) + 1.
            allowed = ip .lt. mr
            p = grofu( real( ip - 1 ) ) / lscale
            j = ( ip - 1 ) * ma + m
            k = j + mr * ma
          else if( lgsp )then
            ip = 4. * p + .1
            allowed = 2 * abs( ip ) .le. mr
            ip = ip + mr / 2 + 1
            j = ( ip - 1 ) * ma + m
            k = j + mr * ma
          else if( s3dc )then
            ip = uofgr( p * lscale ) + 1.
            allowed = ip .le. mr
            p = grofu( real( ip - 1 ) ) / lscale
            j = ( ( m + 1 ) * ( m + 2 ) ) / 2
            j = 4 * ( s3ntm * ( ip - 1 ) + j - 1 ) + 1
            k = j + 1
          else if( scfc )then
            ip = nint( p )
            allowed = ip .le. np
            j = mact * np + 2 * ip - 1
            k = j + 1
          else if( sfpc )then
            ip = nint( p )
            allowed = ip .le. np
            j = 2 * ifun( ip )
            k = j + 1
          else if( sphb )then
            ip = nint( p - 1 )
            allowed = ip .le. ma * ( m + 1 )
            j = ip + ma * ( ( m * ( m + 1 ) ) / 2 )
            k = j + ma * ( ( mr + 1 ) * ( mr + 2 ) ) / 2
            man = ma
          else if( zanl )then
            ip = p * lscale * drfac + .5
            allowed = ip .lt. mr
            p = ( real( ip ) - .5 ) / ( drfac * lscale )
            j = ( ip - 1 ) * ( ma + 1 ) + m + 1
            k = j + mr * ( ma + 1 )
          end if
        end if
c read through file
        if( allowed )then
          do while ( ifail .eq. 0 )
            istp = istp + 1
            if( istp .gt. mstp )call crash( 'LOGANL',
     +                                    'Workspace arrays too small' )
            if( m .gt. 0 )then
              a = sqrt( wres( j ) * wres( j ) + wres( k ) * wres( k ) )
            else
              a = abs( wres( j ) )
            end if
            amp( istp ) = log( max( a, 3.e-7 ) )
            if( ( mact .gt. 0 ) .and. ( a .gt. 0. ) )then
              phase( istp ) =
     +                      atan2( wres( k ), wres( j ) ) / real( mact )
              if( sphb )phase( istp ) = -phase( istp )
            else
              phase( istp ) = 0
            end if
            tm( istp ) = time
            call nextrec( type, ifail )
          end do
c plot growth of selected harmonic
          jstart = max( 1, jstart )
          if( jend .eq. 0 )jend = istp
          if( krun .eq. 1 )then
            call jspage
            call jssize( .1, .7, .1, .9 )
            call jscale( 0., 1., 0., 1. )
            if( nruns .eq. 1 )then
              call jsbldt( 'Run no' )
              call jsbldi( irun, 5 )
            end if
            if( sphb )then
              call jsbldt( ' n=' )
              i = mod( ip - 1, man ) + 1
              call jsbldi( i, 1 )
              call jsbldt( 'm=' )
              i = ( ip - 1 ) / man
              call jsbldi( i, 1 )
              call jsbldt( 'l=' )
              call jsbldi( m, 1 )
            else if( scfc )then
              call jsbldt( ' l=' )
              call jsbldi( mact, 1 )
              man = np / ( mact + 1 )
              i = ( ip - 1 ) / man
              call jsbldt( 'm=' )
              call jsbldi( i, 1 )
              i = ip - 1 - i * man
              call jsbldt( 'n=' )
              call jsbldi( i, 1 )
            else
              call jsbldt( ' m = ' )
              call jsbldi( mact, 2 )
              if( danl .or. zanl )call jsbldt( ' R=' )
              if( s3dc )call jsbldt( ' r=' )
              if( lgsp )call jsbldt( ' p=' )
              if( sfpc )call jsbldt( ' n=' )
              call jsbldf( p, 6, 2 )
            end if
            call jswrit( .2, .8 )
          end if
          call grwplt( mact, istp, amp, phase, tm, jstart, jend )
        end if
      end do
c switch to next run or original if done
      krun = krun + 1
      if( krun .le. nruns )then
        call cmprun( krun )
        go to 5
      else
        call cmprun( 1 )
        if( compt )go to 2
      end if
      go to 1
      end
