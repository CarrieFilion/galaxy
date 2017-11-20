      subroutine contour( lw, lp, m, power )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to contour the data supplied in array x, which results from
c   some form of harmonic analysis of the particle distribution in a
c   simulation.  It can be either the raw data or its power spectrum
c uses JSPLOT graphics routines
      use aarrays
      implicit none
c
c calling arguments
      integer lp, lw, m
      logical power
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
      include 'inc/orbanl.f'
c
c externals
      external akappa, lndbld, trcont
      character*4 datnam
      logical gtlogl
      real grofu
c
c local arrays
      real c( 10 )
      real, allocatable :: x( :, : )
c
c local variables
      integer i, ip, ires, iw, jlw, llw, n, nc
      logical numer
      real am, prob, r, t, xl, xmin, xmax, xu, y, ymin, ymax
      real x1, x2, y1, y2
      include 'inc/pi.f'
c
      if( .not. ( danl .or. dnst .or. lgsp .or. s3dc .or. sfpc .or.
     +             scfc .or.sphb .or. vfld .or. zanl )
     +                )call crash( 'CONTOUR', 'Unrecognised data type' )
      if( power )then
        print *, 'Entering 0 in answer to this question will include' //
     +           ' negative frequencies'
        print *, '1 will show all frequencies as positive'
        print *, 'a number >1 will show that fraction of the range'
      end if
      call gtintg( 'Enter fraction of array to be contoured', llw )
      jlw = 1
      if( llw .gt. 0 )then
        llw = 1 + ( lw - 1 ) / llw
      else
        llw = lw
        if( power )jlw = -llw / 3
      end if
      call jspage
      numer = .false.
c set scaling for contour routine
      if( danl .or. dnst .or. s3dc .or. sfpc )then
        xl = jp - 1
        xu = kp - 1
      else if( lgsp )then
        xl = .25 * real( jp - ncp / 2 - 1 )
        xu = .25 * real( kp - ncp / 2 - 1 )
      else if( scfc .or. sphb )then
        xl = jp
        xu = kp
      else if( vfld .or. zanl )then
        xl = ( real( jp ) - .5 ) / ( drfac * lscale )
        xu = ( real( kp ) - .5 ) / ( drfac * lscale )
      end if
      if( power )then
        ymin = 2. * pi * real( jlw - 1 ) / ( tme( kt ) - tme( jt ) )
        ymax =
     +      2. * pi * real( llw - 1 ) / ( tme( kt ) - tme( jt ) ) + ymin
        if( danl .or. vfld )numer = .not. gtlogl(
     +                         'Use theoretical frequency curves?' )
      else
        ymin = tme( jt )
        ymax = tme( kt )
      end if
      xscale = ( xu - xl ) / real( lp - 1 )
      yscale = ( ymax - ymin ) / real( llw - 1 )
      xoffs = xl - xscale
      yoffs = ymin - yscale
c set up frame
      if( danl .or. dnst .or. s3dc .or. vfld .or. zanl )then
        xmin = 0
        xmax = rgrid( 1 ) / lscale
      else if( lgsp )then
        xmax = .25 * real( ncp / 2 )
        xmin = -xmax
      else if( sfpc )then
        xmin = 0
        xmax = ncp - 1
      else if( scfc .or. sphb )then
        xmin = 1
        xmax = ncp
      end if
      call jssize( .1, .9, .1, .8 )
      call jscale( xmin, xmax, ymin, ymax )
      call jsthik( 3 )
c draw radial limits of data used
      if( danl .or. dnst .or. s3dc )then
        xl = grofu( xl ) / lscale
        xu = grofu( xu ) / lscale
      end if
      call jsdash( 0, 2, 0, 2 )
      if( xl .gt. xmin )then
        call jsmove( xl, ymin )
        call jsline( xl, ymax )
      end if
      if( xu .lt. xmax )then
        call jsmove( xu, ymin )
        call jsline( xu, ymax )
      end if
      call jsdash( 0, 0, 0, 0 )
c write heading
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      if( sphb )then
        call jsbldt( ' l =' )
      else
        call jsbldt( ' m =' )
      end if
      call jsbldi( m, 2 )
      call jswrit( xmin, 1.15 * ymax - .15 * ymin )
c assemble information for sub-heading
      do i = 1, ndattyp
        if( dattyp( i ) )n = i
      end do
      call jsbldt( datnam( n ) )
      call jsbldt( 'data' )
      if( power )then
        call jsbldt( '  over time range' )
        call jsbldf( tme( jt ), 5, 1 )
        call jsbldt( 'to' )
        call jsbldf( tme( kt ), 5, 1 )
      end if
      if( apodc .ne. 0. )then
        call jsbldt( '  scaling exponent' )
        call jsbldf( apodc, 6, 3 )
      end if
      call jschsz( .2 )
      call jswrit( .9 * xmin + .1 * xmax, 1.08 * ymax - .08 * ymin )
      call jschsz( .3 )
c draw and label axes
      if( lgsp )call jsaxis( 'X', 'tan(\\gamma)', 1 )
      if( danl .or. dnst .or. s3dc .or. vfld .or. zanl )then
        if( scale_set )then
          x1 = xmin * unit_L
          x2 = xmax * unit_L
          call jscale( x1, x2, ymin, ymax )
          call jsaxis( 'X', 'R (kpc)', 1 )
          call jscale( xmin, xmax, ymin, ymax )
        else
          call jsaxis( 'X', 'radius', 1 )
        end if
      end if
      if( sfpc .or. scfc )call jsaxis( 'X', 'Basis function', 1 )
      if( sphb )call jsaxis( 'X', 'Spherical Bessel function', 1 )
      if( power )then
        if( scale_set )then
          y1 = ymin * unit_V / unit_L
          y2 = ymax * unit_V / unit_L
          call jscale( xmin, xmax, y1, y2 )
          call jsaxis( 'y', 'km s^{-1} kpc^{-1}', 1 )
          call jscale( xmin, xmax, ymin, ymax )
        else
          call jsaxis( 'y', 'Frequency', 1 )
        end if
      else
        if( scale_set )then
          y1 = ymin * unit_T
          y2 = ymax * unit_T
          call jscale( xmin, xmax, y1, y2 )
          call jsaxis( 'y', 'Time (Myr)', 1 )
          call jscale( xmin, xmax, ymin, ymax )
        else
          call jsaxis( 'Y', 'Time', 1 )
        end if
      end if
c pick out data
      allocate ( x( lp, lw ) )
      n = 0
      am = 0.
      if( .not. power )llw = lw
      do iw = 1, llw
        do ip = 1, lp
          if( ( .not. power ) .and. ( m .eq. 0 ) )then
            n = n + 1
            x( ip, iw ) = sdata( 1, n )
          else
            nc = iw + jlw - 1
            if( nc .gt. 0 )then
              n = n + 1
              x( ip, iw ) = sdata( 1, n )**2 + sdata( 2, n )**2
            else
              nc = ip + lp * ( lw + nc - 1 )
              x( ip, iw ) = sdata( 1, nc )**2 + sdata( 2, nc )**2
            end if
          end if
          if( iw .gt. 1 )am = max( am, x( ip, iw ) )
        end do
      end do
c choose contour levels
      nc = 10
      if( .not. power )nc = 5
      prob = 0
c      prob = pi / real( 4 * nbod )
c      prob = sqrt( prob )
      if( power .or. ( am .lt. prob ) )then
        do i = 1, nc
          c( i ) = .5 * am * ( real( i ) - .5 ) / real( nc )
        end do
      else
        do i = 1, nc
          c( i ) = .5 * am * ( real( i ) - .5 ) / real( nc )
c          c( i ) = prob * real( 2 * i )
c          c( i ) = .04 * real( i )
        end do
c        print *, 'contour levels reset to .04*i'
      end if
c contour array
      call jsthik( 1 )
      call jscont( x, lp, llw, c, nc, trcont )
      call jsthik( 3 )
c
c      prob = ( tme( kt ) - tme( jt ) ) / ( 2. * pi )
c      call gtreal( 'Enter frequency', am )
c      do while ( am .gt. 0. )
c        i = am * prob
c        i = max( 0, i )
c        i = min( i, nt - 1 )
c        am = real( i ) / prob
c        call jsmove( xmin, am )
c        call jsline( xmax, am )
c        call gtreal( 'Enter frequency', am )
c      end do
c
      if( power )then
        if( danl .or. vfld )then
c draw grid lines
c          call jsdash( 0, 3, 0, 3 )
c          do i = 1, nr
c            xx = i - 1
c            xx = grofu( xx ) / lscale
c            if( xx .gt. xmin .and. xx .lt. xmax )then
c              call jsmove( xx, 0. )
c              call jsline( xx, ymax )
c            end if
c          end do
c draw resonance lines
          if( numer )then
c read frequency data
            call rewnd
            call nextrec( 'FRQS', i )
            do while ( i .eq. 0 )
              t = ts * real( istep )
              t = min( abs( t - tme( jt ) ), abs( t - tme( kt ) ) )
              if( t .lt. .1 )then
                if( danl )then
c work over resonances
                  do ires = 1, 3
                    if( ires .eq. 2 )then
                      call jsdash( 0, 0, 0, 0 )
                    else
                      call jsdash( 2, 2, 2, 2 )
                    end if
c work over radii
                    do i = 2, mr - 1
                      r = ( rinh( 1 ) + real( i ) * drfrqs ) / lscale
                      y = real( m ) * wres( i )
                      if( ires .eq. 1 )y = y + wres( i + mr )
                      if( ires .eq. 3 )y = y - wres( i + mr )
                      y = min( y, ymax )
                      if( i .eq. 2 )call jsmove( r, y )
                      call jsline( r, y )
                    end do
                  end do
                else
                  call jsdash( 2, 2, 2, 2 )
c work over radii
                  do i = 2, mr - 1
                    r = ( rinh( 1 ) + real( i ) * drfrqs ) / lscale
                    y = wres( i + mr )
                    y = min( y, ymax )
                    if( i .eq. 2 )call jsmove( r, y )
                    call jsline( r, y )
                  end do
                end if
              end if
              call nextrec( 'FRQS', i )
            end do
            call jsdash( 0, 0, 0, 0 )
          else
c use theoretical frequencies
            if( danl )then
              morb = m
              do ires = 1, 3
                lorb = ires - 2
                if( ires .eq. 2 )then
                  call jsdash( 0, 0, 0, 0 )
                else
                  call jsdash( 2, 2, 2, 2 )
                end if
                call jsplt2( lndbld )
              end do
            else
              call jsdash( 2, 2, 2, 2 )
              call jsplt2( akappa )
            end if
          end if
          call jsdash( 0, 0, 0, 0 )
        end if
      end if
c finish up
      call jschsz( .3 )
      return
      end
