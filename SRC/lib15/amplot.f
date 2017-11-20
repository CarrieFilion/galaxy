      subroutine amplot
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot results from several different types of analysis saved
c   in the runXXX.res file.  The variation of the mean of a range of values
c   as a function of time is plotted
c
c Data types handled by this routine are:
c      danl - sectoral harmonics of the mass distbn on grid rings
c      lgsp - logarithmic spiral coeffs
c      sfpc - coeffs saved from the SFP expansion
c      sphb - spherical Bessel functions coeffs
c      zanl - sectoral harmonics of the mid-plane displacements
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
      include 'inc/jscmmn.f'
c
      include 'inc/model.f'
c
c externals
      logical gtlogl
      real grofu, roundup, uofgr
c      real*8 gsigma
c
c local arrays
      integer ifun( mostn + 1 )
c      integer isize
c      parameter ( isize = 5001 )
      real, allocatable :: amp( : ), tm( : )
      real aa( 2 ), bb( 2 )
c
c local variables
      character prompt*80, type*4
      integer i, ifail, ip, isize, istp, j, jcmp, k, krun, l, m, mact
      logical logy, ok, phase, skip
      real a, area, b, c, p, r, r1, r2, xmax, ymax, ymin
c      real*8 r2
      include 'inc/pi.f'
c
c select type of data
      call datype( type )
      jcmp = icmp
c count number of available records and allocate space
      call rewnd
      call nextrec( type, ifail )
      isize = 0
      do while ( ifail .eq. 0 )
        isize = isize + 1
        call nextrec( type, ifail )
      end do
      call rewnd
      allocate ( amp( isize ) )
      allocate ( tm( isize ) )
c
      logy = gtlogl( 'Logarithmic amplitude scale?' )
c choose sectoral harmonic
      call choosm( type, mact, m )
      do while ( mact .ge. 0 )
        if( sfpc )then
          ncp = 0
          do i = 1, lastf
            if( mact .eq. msel( i ) )then
              ncp = ncp + 1
              if( ncp .gt. mostn + 1 )call crash( 'AMPLOT',
     +                                          'ifun array too small' )
              ifun( ncp ) = i
            end if
          end do
          if( ncp .eq. 0 )call crash( 'AMPLOT',
     +                        'No basis function with this m selected' )
        end if
c
        phase = .false.
        if( ( mact .gt. 0 ) .and. ( danl .or. zanl .or. sphb ) )
     +                     phase = gtlogl( 'Discard phase variations?' )
c select radius or whatever to plot
        ok = .false.
        do while ( .not. ok )
          if( danl .or. zanl )write( prompt, '( a )' )
     +                              'Select inner radius or -ve to skip'
          if( lgsp )write( prompt, '( a )' )
     +                  'Select tan(gamma) value or large value to skip'
          if( s3dc )write( prompt, '( a )' )
     +                                    'Select radius or -ve to skip'
          if( sfpc )then
            j = ncp - 1 + nsel( ifun( 1 ) )
            write( prompt, '( a, i3, a)' )
     +      'Select first basis function (up to', j, '), or -ve to skip'
          end if
          if( sphb )write( prompt, '( a )' )
     +                                  'Select function or -ve to skip'
          call gtreal( prompt, p )
          if( lgsp )then
            skip = abs( p ) .gt. 100.
          else
            skip = p .lt. 0.
          end if
          if( .not. skip )then
c identify required data locations
            call rewnd
            call nextrec( type, ifail )
            if( danl )then
              jp = uofgr( p * lscale ) + 1.5
              ok = jp .lt. mr
              if( ok )then
                jp = max( jp, 1 )
                p = grofu( real( jp - 1 ) ) / lscale
                call gtreal(
     +              'Select outer radius - values will be averaged', r )
                kp = uofgr( r * lscale ) + 1.5
                r = grofu( real( kp - 1 ) ) / lscale
                print *, 'Radial range', p, ' to', r, ' selected'
                ok = ( kp .lt. mr ) .and. ( kp .ge. jp )
              end if
            else if( lgsp )then
              ip = 4. * p + .1
              ok = 2 * abs( ip ) .le. mr
              if( ok )then
                ip = ip + mr / 2 + 1
                j = ( ip - 1 ) * ma + m
                k = j + mr * ma
              end if
            else if( s3dc )then
              jp = uofgr( p * lscale ) + 1.5
              ok = jp .lt. mr
              if( ok )then
                jp = max( jp, 1 )
                p = grofu( real( jp - 1 ) ) / lscale
                j = ( ( m + 1 ) * ( m + 2 ) ) / 2
                j = 4 * ( s3ntm * ( jp - 1 ) + j - 1 ) + 1
                k = j + 1
              end if
            else if( sfpc )then
              ip = nint( p ) + 1 - nsel( ifun( 1 ) )
              ok = ( ip .gt. 0 ) .and. ( ip .le. ncp )
            else if( sphb )then
              jp = nint( p )
              k = ma * ( m + 1 )
              ok = jp .le. k
              if( ok )then
                call gtreal(
     +             'Select last function - values will be averaged', a )
                kp = nint( a )
                print *, 'Functions', jp, ' to', kp, ' selected'
                ok = kp .le. k
                if( ( jp - 1 ) / ma .ne. ( kp - 1 ) / ma )then
                  print *, 'Odd/even functions mixed'
                  ok = .false.
                end if
              end if
            else if( zanl )then
              jp = p * lscale * drfac + .5
              ok = jp .lt. mr
              if( ok )then
                jp = max( jp, 1 )
                p = ( real( jp ) - .5 ) / ( drfac * lscale )
                call gtreal(
     +             'Select outer radius - values will be averaged', r )
                kp = r * lscale * drfac + .5
                r = ( real( kp ) - .5 ) / ( drfac * lscale )
                print *, 'Radial range', p, ' to', r, ' selected'
                ok = ( kp .lt. mr ) .and. ( kp .ge. jp )
              else
                call crash( 'AMPLOT', 'Unrecognized data type' )
              end if
            end if
c label plot
            if( ok )then
              call jspage
              call jssize( 0., 1., 0., 1. )
              call jscale( 0., 1., 0., 1. )
              if( nruns .eq. 1 )then
                call jsbldt( 'Run no' )
                call jsbldi( irun, 5 )
              end if
              call jsbldt( ' l =' )
              call jsbldi( mact, 2 )
              if( danl )call jsbldt( ' density amplitude' )
              if( lgsp )call jsbldt( ' logarithmic spiral' )
              if( s3dc )call jsbldt( ' surface harmonic' )
              if( sfpc )call jsbldt( ' basis function' )
              if( sphb )call jsbldt( ' spherical Bessel' )
              if( zanl )call jsbldt( ' bending amplitude' )
              call jswrit( .2, .9 )
c
              if( danl .or. zanl )then
                call jsbldf( p, 5, 2 )
                call jsbldt( '< R <' )
                call jsbldf( r, 5, 2 )
              else if( lgsp )then
                call jsbldt( 'tan\\gamma=' )
                call jsbldf( p, 5, 2 )
              else if( s3dc )then
                call jsbldt( 'r =' )
                call jsbldf( p, 5, 2 )
              else if( sfpc )then
                call jsbldt( 'n=' )
                call jsbldi( ip - 1, 2 )
              else if( sphb )then
                call jsbldt( 'm=' )
                i = ( jp - 1 ) / ma
                call jsbldi( i, 1 )
                call jsbldt( 'n=' )
                i = mod( jp - 1, ma ) + 1
                call jsbldi( i, 1 )
                if( kp .gt. jp )then
                  call jsbldt( 'to' )
                  i = ( kp - 1 ) / ma
                  if( i .gt. ( jp - 1 ) / ma )then
                    call jsbldt( 'm=' )
                    call jsbldi( i, 1 )
                    call jsbldt( 'n=' )
                  end if
                  i = mod( kp - 1, ma ) + 1
                  call jsbldi( i, 1 )
                end if
              end if
              call jswrit( .2, .75 )
c
              do krun = 1, nruns
                istp = 0
                if( krun .eq. 1 )then
c select grid appropriate for population
                  icmp = jcmp
                  if( hybrid )then
                    i = igrd( icmp )
                    call switch( i )
                  end if
                  ymin = 0
                  ymax = 0
                end if
c work through file
                do while ( ifail .eq. 0 )
                  istp = istp + 1
                  tm( istp ) = time
                  a = 0
                  b = 0
                  c = 0
                  if( danl .or. sphb .or. zanl )then
                    if( danl )then
                      r1 = max( real( jp ) - .5, 0. )
                      r1 = grofu( r1 ) / lscale
                    end if
                    do ip = jp, kp
                      if( danl )then
c average weighted by mass not surface density
                        r2 = ip
                        r2 = grofu( r2 - .5 ) / lscale
c .5 * alpha = pi / na
                        area = .5 * alpha * ( r2 * r2 - r1 * r1 )
                        r1 = r2
                        l = ( ip - 1 ) * ma + 1
                        wres( l ) = wres( l ) * area
                        j = l + m
                        k = j + mr * ma
                        if( m .gt. 0 )then
                          c = c + wres( l )
                          wres( j ) = wres( l ) * wres( j )
                          wres( k ) = wres( l ) * wres( k )
                        end if
                      else if( s3dc )then
                        j = ( ip - 1 ) * ( ma + 1 ) + m + 1
                        k = j + mr * ( ma + 1 )
                      else if( sphb )then
                        j = ip + ma * ( ( m * ( m + 1 ) ) / 2 )
                        k = j + ma * ( ( mr + 1 ) * ( mr + 2 ) ) / 2
                      else if( zanl )then
                        j = ( ip - 1 ) * ( ma + 1 ) + m + 1
                        k = j + mr * ( ma + 1 )
                      end if
                      if( mact .eq. 0 )then
                        a = a + wres( j )
                      else
                        if( phase )then
                          a = a + wres( j )
                          b = b + wres( k )
                        else
                          a = a + sqrt( wres( j ) * wres( j ) +
     +                         wres( k ) * wres( k ) )
                        end if
c                        if( danl )then
c                          r2 = p
c                          a = a * wres( j - m + 1 ) /
c     +                              ( real( na ) * gsigma( r2 ) )
c                        end if
                      end if
                    end do
                    if( phase )a = sqrt( a * a + b * b )
                    if( danl )then
                      if( ( m .gt. 0 ) .and. ( c .gt. 0. ) )a = a / c
                    else
                      a = a / real( kp - jp + 1 )
                    end if
                    amp( istp ) = a
                  else
                    if( sfpc )then
                      j = 2 * ifun( ip )
                      k = j + 1
                    end if
                    if( mact .eq. 0 )then
                      amp( istp ) = wres( j )
                    else
                      a = sqrt( wres( j ) * wres( j ) +
     +                          wres( k ) * wres( k ) )
c
c                    if( lgsp .and. uqmass .and.
c     +                  ( irun .lt. 465 ) )a = a / fmass( icmp )
c
                      amp( istp ) = a
                    end if
                  end if
                  if( logy )then
                    a = max( abs( amp( istp ) ), 1.e-10 )
                    amp( istp ) = log( a )
                  end if
                  if( krun .eq. 1 )then
                    ymin = min( ymin, amp( istp ) )
                    ymax = max( ymax, amp( istp ) )
                  end if
                  call nextrec( type, ifail )
                  if( ( ifail .eq. 0 ) .and. ( istp .eq. isize ) )then
                    print *, 'Local arrays too small in AMPLOT'
                    ifail = 1
                  end if
                end do
c draw axes
                if( krun .eq. 1 )then
                  xmax = ( int( .1 * tm( istp ) ) + 1 ) * 10
                  ymax = roundup( ymax )
                  ymin = -roundup( -ymin )
                  call jssize( .15, .95, .15, .85 )
                  if( scale_set )then
                    call jscale( 0., xmax * unit_T, ymin, ymax )
                    call jsaxis( 'X', 'Time (Myr)', 4 )
                    call jscale( 0., xmax, ymin, ymax )
                  else
                    call jscale( 0., xmax, ymin, ymax )
                    call jsaxis( 'X', 'Time', 4 )
                  end if
                  if( logy )then
                    call jsaxis( 'Y', 'ln(A)', 1 )
                  else
                    call jsaxis( 'Y', 'A', 1 )
                  end if
                  if( lgsp )then
                    a = pi / ( 4. * nbod )
                    a = sqrt( a )
                    if( logy )a = log( a )
                    call jsdash( 2, 2, 2, 2 )
                    call jsmove( 0., a )
                    call jsline( xmax, a )
                    call jsdash( 0, 0, 0, 0 )
                  end if
                end if
c plot data
                call jscoin( tm, amp, istp, krun )
c mark color code
                if( nruns .gt. 1 )then
                  aa( 1 ) = .05 * xmax
                  b = .04 * real( krun )
                  bb( 1 ) = b * ymax + ( 1. - b ) * ymin
                  aa( 2 ) = aa( 1 ) + .05 * xmax
                  bb( 2 ) = bb( 1 )
                  call jscoin( aa, bb, 2, krun )
                  call jsbldi( irun, 5 )
                  a = aa( 2 ) + .01 * xmax
                  b = bb( 2 ) - .01 * ( ymax - ymin )
                  call jswrit( a, b )
                end if
c switch to next run or original if done
                if( nruns .gt. 1 .and. krun .lt. nruns )then
                  call cmprun( krun + 1 )
                  call rewnd
                  call nextrec( type, ifail )
c allow for different nsect, mr, and ma
                  m = mact / nsect
                  if( lgsp )then
                    ip = 4. * p + .1
                    ip = ip + mr / 2 + 1
                    j = ( ip - 1 ) * ma + m
                    k = j + mr * ma
                  end if
                end if
c end krun loop
              end do
              if( nruns .gt. 1 )call cmprun( 1 )
            end if
c end skip block
          end if
c skip = T flags end of this sectoral harmonic
          ok = skip
c end ok loop
        end do
c input new sectoral harmonic
        call choosm( type, mact, m )
      end do
      return
      end
