      subroutine slowng
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot results from several different types of analysis saved
c   in the runXXX.res file.  The variation at specific times is plotted
c
c Data types handled by this routine are:
c      danl - sectorial harmonics of the mass distbn on grid rings
c      dnst - radial mass distribution on grid rings
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
      include 'inc/grids.f'
c
c externals
      real grofu, roundup, uofgr
c
c local arrays
      integer ifun( mostn + 1 )
      integer mstp
      parameter (mstp = 5001 )
      real amp( mstp ), phase( mstp )
c
c local variables
      character type*4
      integer i, ifail, ip, istp, j, k, m, mact, man, np
      logical allowed
      real a, p, pif, pit, xmax, ymax, ymin
      include 'inc/pi.f'
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
            if( np .gt. mostn + 1 )call crash( 'SLOWNG',
     +                                          'ifun array too small' )
            ifun( np ) = ip
          end if
        end do
        if( np .eq. 0 )call crash( 'SLOWNG',
     +                        'No basis function with this m selected' )
      end if
c select p value for analysis
      allowed = mact .ne. 0
      do while ( allowed )
c rewind file
        call rewnd
        call nextrec( type, ifail )
        if( .not. ( scfc .or. sfpc ) )np = mr
        if( danl .or. zanl )call gtreal( 'Select radius', p )
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
c set pointers to pick out requested values
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
c rewind and read through file
        call rewnd
        call nextrec( type, ifail )
        istp = 0
        if( allowed )then
          do while ( ifail .eq. 0 )
            istp = istp + 1
            if( istp .gt. mstp )call crash( 'SLOWNG',
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
            else
              phase( istp ) = 0
            end if
            tme( istp ) = time
            call nextrec( type, ifail )
          end do
c check for sufficient data
          if( istp .le. 16 )then
            print *, istp, ' data values available, need at least 17'
            return
          end if
c make phases monotonic
          pif = 0.
          pit = 0.
          if( mact .ne. 0 )pit = 2. * pi / real( mact )
          do i = 2, istp
            if( ( phase( i ) + pif ) .lt. ( phase( i - 1 ) - .5 * pit )
     +          )pif = pif + pit
            phase( i ) = phase( i ) + pif
          end do
c get slope of phase at each time - Lagrange formula for uniformly spaced data
          a = ( tme( istp ) - tme( 1 ) ) / real( istp - 1 )
          do i = 5, istp - 4
c eight-point formula for first derivative
            phase( i - 4 ) =
     +     ( 672. * ( phase( i + 1 ) - phase( i - 1 ) ) -
     +       168. * ( phase( i + 2 ) - phase( i - 2 ) ) +
     +        32. * ( phase( i + 3 ) - phase( i - 3 ) ) -
     +         3. * ( phase( i + 4 ) - phase( i - 4 ) ) ) / ( 840. * a )
          end do
c smoothing
          do i = 1, istp - 16
            a = 0
            do k = 1, 9
              a = a + phase( i + k - 1 )
            end do
            phase( i ) = a / 9.
          end do
c set window
          call jspage
          call jssize( .15, .95, .15, .95 )
          xmax = roundup( tme( istp ) )
          ymin = 0
          ymax = 0
          do i = 1, istp - 16
            ymin = min( ymin, phase( i ) )
            ymax = max( ymax, phase( i ) )
          end do
          ymax = roundup( ymax )
          call jscale( 0., xmax, ymin, ymax )
c draw axes and plot pattern speed of selected harmonic
          call jsaxis( 'x', 'time', 1 )
          call jsaxis( 'y', '\\Omega_{p}', 1 )
c label plot
          call jsbldt( 'Run no' )
          call jsbldi( irun, 5 )
          if( sphb )then
            call jsbldt( ' n =' )
            i = mod( ip - 1, man ) + 1
            call jsbldi( i, 1 )
            call jsbldt( 'm =' )
            i = ( ip - 1 ) / man
            call jsbldi( i, 1 )
            call jsbldt( 'l =' )
            call jsbldi( m, 1 )
          else if( scfc )then
            call jsbldt( ' l =' )
            call jsbldi( mact, 1 )
            man = np / ( mact + 1 )
            i = ( ip - 1 ) / man
            call jsbldt( 'm =' )
            call jsbldi( i, 1 )
            i = ip - 1 - i * man
            call jsbldt( 'n =' )
            call jsbldi( i, 1 )
          else
            call jsbldt( ' m =' )
            call jsbldi( mact, 2 )
            if( danl .or. zanl )call jsbldt( ' r =' )
            if( lgsp )call jsbldt( ' tan\\gamma =' )
            if( sfpc )call jsbldt( ' n=' )
            call jsbldf( p, 6, 2 )
          end if
          call jswrit( .2 * xmax, 1.05 * ymax - 0.05 * ymin )
          end if
          call jsjoin( tme( 9 ), phase, istp - 16 )
      end do
      go to 1
      end
