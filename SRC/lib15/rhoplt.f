      subroutine rhoplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot radial density profiles saved in the runXXX.res file.
c    The variation at specific times is plotted
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
      external gmassd, gmassh, lgmassh, lgsigma, lrhohal
      real roundup
      real*8 gmassd, gmassh
c
c local variables
      character string*4, type*4
      integer i, iend, ifail, istp, j, jstep, krun
      logical first, logm, rhop
      real mfac, r1, r2, t, tmax, x, xlast, xmax, xmin, xstep, y, ymax
      real ymin, z, x1, x2, y1, y2
      include 'inc/pi.f'
c      data logm / .true. /
c
c determine which population to plot
      call selpop( icmp )
      i = igrd( icmp )
      call switch( i )
      if( disc( icmp ) )then
        type = 'SIGR'
        logm = .false.
      else
        type = 'RHOR'
        logm = .true.
      end if
c
      call rewnd
      call nextrec( type, ifail )
      if( ifail .ne. 0 )then
        print *, 'No data found for type ', type
        return
      end if
c
    1 call gtchar( 'Plot density or mass profile?', string )
      call uppercase( string )
      rhop = string .eq. 'DENS'
      if( ( .not. rhop ) .and. ( string .ne. 'MASS' ) )go to 1
c choose plotting times
      call lstrid( type, jstep, iend )
c choose scaling
      if( rhop )then
        if( disc( icmp ) )then
          xmin = 0
          xmax = max( rgrid( jgrid ) / lscale, wres( mr - 1 ) )
          xmax = nint( xmax )
          ymin = -10
          if( wres( mr ) .gt. 0. )ymin = log( wres( mr ) ) + 1.
          ymax = 0
          if( scale_set )then
            x1 = xmin * unit_L
            x2 = xmax * unit_L
            y1 = log( 1.e+4 * exp( ymin ) * unit_M / unit_L**2 )
            y2 = log( 1.e+4 * exp( ymax ) * unit_M / unit_L**2 )
          end if
        else
          xmin = log10( wres( 1 ) ) - 1.
          xmax = log10( wres( mr - 1 ) ) + 1.
          xmax = nint( xmax )
          xmin = -roundup( -xmin ) + 1
          ymin = log10( wres( mr ) )
          ymin = nint( ymin ) - 1
          ymax = log10( wres( 2 ) ) + 2.
          ymax = roundup( ymax )
          if( scale_set )then
            x1 = log10( 10.**xmin * unit_L )
            x2 = log10( 10.**xmax * unit_L )
            y1 = log10( 10.**ymin * unit_M / unit_L**3 ) + 1.
            y2 = log10( 10.* ymax * unit_M / unit_L**3 ) + 1.
          end if
        end if
      else
        xmin = 0
        xmax = max( rgrid( jgrid ) / lscale, wres( mr - 1 ) )
        ymin = 0
        if( ctype( icmp ) .ne. 'UNKN' )then
          if( disc( icmp ) )then
            ymax = gmassd( dble( xmax ) )
            if( testp( icmp ) )ymax = 1
          else
            ymax = gmassh( dble( xmax ) )
          end if
c          ymax = roundup( .99 * ymax )
        else
          y = 0
          r2 = 0
          do i = 1, mr, 2
            r1 = r2
            r2 = 2 * wres( i ) - r1
            if( disc( icmp ) )then
              y = y + pi * wres( i + 1 ) * ( r2**2 - r1**2 )
            else
              y = y + 4. * pi * wres( i + 1 ) * ( r2**3 - r1**3 ) / 3.
            end if
          end do
          ymax = roundup( y )
        end if
        if( logm )then
          xmax = log10( xmax )
          xmin = xmax - 5
          ymax = log10( ymax )
          ymin = ymax - 5
          if( scale_set )then
            x1 = log10( 10.**xmin * unit_L )
            x2 = log10( 10.**xmax * unit_L )
            y1 = log10( 10.**ymin * unit_M )
            y2 = log10( 10.* ymax * unit_M )
          end if
        else
          if( scale_set )then
            x1 = xmin * unit_L
            x2 = xmax * unit_L
            y1 = ymin * unit_M
            y2 = ymax * unit_M
          end if
        end if
      end if
c draw and mark axes
      call jspage
      call jssize( .15, .7, .1, .85 )
      if( scale_set )then
        call jscale( x1, x2, y1, y2 )
      else
        call jscale( xmin, xmax, ymin, ymax )
      end if
      if( rhop )then
        if( disc( icmp ) )then
          if( scale_set )then
            call jsaxis( 'x', 'R (kpc)', 1 )
            call jsaxis( 'y', 'ln\\Sigma (M_0 pc^{-2})', 1 )
          else
            call jsaxis( 'x', 'R', 1 )
            call jsaxis( 'y', 'ln\\Sigma', 1 )
          end if
        else
          if( scale_set )then
            call jslaxs( 'x', 'r (kpc)', 1 )
            call jslaxs( 'y', '\\rho (M_0 pc^{-3})', 1 )
          else
            call jslaxs( 'x', 'r', 1 )
            call jslaxs( 'y', '\\rho', 1 )
            if( s3d )then
              x = log10( .5 * s3rad( 2 ) )
              if( x .gt. xmin )then
                call jsdash( 0, 1, 0, 1 )
                call jsmove( x, ymin )
                call jsline( x, ymax )
                call jsdash( 0, 0, 0, 0 )
              end if
            end if
          end if
        end if
      else
        if( disc( icmp ) )then
          if( scale_set )then
            call jsaxis( 'x', 'R (kpc)', 1 )
            call jsaxis( 'y', 'M(R)/10^{10}', 1 )
          else
            call jsaxis( 'x', 'R', 1 )
            call jsaxis( 'y', 'M(R)', 1 )
          end if
        else
          if( logm )then
            call jslaxs( 'x', 'r', 1 )
            call jslaxs( 'y', 'M(r)', 1 )
          else
            call jsaxis( 'x', 'r', 1 )
            call jsaxis( 'y', 'M(r)', 1 )
          end if
        end if
      end if
      if( scale_set )call jscale( xmin, xmax, ymin, ymax )
c work through file
      krun = 1
    5 if( nruns .gt. 1 )call pgsci( krun + 1 )
      ifail = 0
      istp = 0
      do while ( ifail .eq. 0 )
        istp = istp + 1
c plot density variation at this moment - reconstructing geometric mean radius
        if( rhop )then
c reconstruct largest particle radius and mass in each bin
          r2 = 0
          do i = 1, mr, 2
            r1 = r2
            r2 = 2 * wres( i ) - r1
            wres( i ) = r2
            if( disc( icmp ) )then
              x = pi * ( r2**2 - r1**2 )
            else
              x = 4. * pi * ( r2**3 - r1**3 ) / 3.
            end if
            wres( i + 1 ) = x * wres( i + 1 )
          end do
c combine bins to improve S/N while retaining resolution
          r2 = 0
          i = -1
          first = .true.
          do while ( i + 2 .lt. mr )
            r1 = r2
            i = i + 2
            r2 = wres( i )
            y = wres( i + 1 )
            j = 1
            do while ( ( r2 .lt. 1.05 * r1 ) .and. ( i + 2 .lt. mr ) )
c     +              .or. ( j .le. 20 ) )
              j = j + 1
              if( i + 2 .lt. mr )then
                i = i + 2
                r2 = wres( i )
                y = y + wres( i + 1 )
              end if
            end do
            if( y .gt. 0. )then
              if( disc( icmp ) )then
c reconstruct surface density
                x = pi * ( r2**2 - r1**2 )
                y = log( y / x )
c reconstruct mean radius
                x = sqrt( .5 * ( r2**2 + r1**2 ) )
              else
c reconstruct volume density
                x = 4. * pi * ( r2**3 - r1**3 ) / 3.
                y = y / x
                y = log10( y )
c reconstruct mean radius
                x = .5 * ( r2**3 + r1**3 )
                x = log10( x ) / 3.
              end if
              if( first )call jsmove( x, y )
              first = .false.
              call jsline( x, y )
            end if
          end do
        else
          r2 = 0
          if( disc( icmp ) )then
            mfac = fmass( 1 )
          else
            mfac = fmass( 2 )
          end if
c          mfac = mfac * ( dmass + hmass ) * real( ma ) / real( nbod )
          mfac = mfac * cmpmas( icmp ) * real( ma ) / real( nbod )
          y = 0
          x = 0
          xstep = .01 * ( xmax - xmin )
          xlast = x
          do i = 1, mr, 2
            r1 = r2
            r2 = 2 * wres( i ) - r1
            x = r2
            if( disc( icmp ) )then
              y = y + pi * wres( i + 1 ) * ( r2**2 - r1**2 )
              if( logm )then
                x = log10( x )
                z = log10( y )
              else
                z = y
              end if
              if( i .eq. 1 )call jsmove( x, z )
              call jsline( x, z )
            else
              y = y + 4. * pi * wres( i + 1 ) * ( r2**3 - r1**3 ) / 3.
              if( logm )then
                x = log10( x )
                z = log10( y )
              else
                z = y
              end if
              if( x .lt. -2.4 )xlast = x
              if( i .eq. 1 .or. x .eq. xlast )call jsmove( x, z )
c              if( i .eq. 1 )call jsmove( x, z )
              if( abs( x - xlast ) .gt. xstep )then
                call jsline( x, z )
                xlast = x
              end if
            end if
c            y = .5 * mfac * real( i + 1 )
c            if( i .eq. 1 )call jsmove( x, y )
          end do
        end if
c get next record needed record
        call nextrec( type, ifail )
        if( istep .gt. iend )ifail = 1
        do while ( ( ifail .eq. 0 ) .and.
     +             ( mod( istep, jstep ) .ne. 0 ) )
          call nextrec( type, ifail )
          if( istep .gt. iend )ifail = 1
        end do
        if( krun .eq. 1 )then
          tmax = time
        else
          if( time .gt. tmax )ifail = 1
        end if
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
        call jsbldt( ' time interval' )
        x = ts
        i = jstep * ( istep / jstep )
        if( istp .gt. 1 )x = ts * real( i ) / real( istp - 1 )
        call jsbldf( x, 6, 1 )
        x = .8 * xmin + .2 * xmax
        y = -.05 * ymin + 1.05 * ymax
        call jswrit( x, y )
c draw theoretical curve
        if( ctype( icmp ) .ne. 'UNKN' )then
          call jsdash( 0, 1, 0, 1 )
          if( rhop )then
            if( disc( icmp ) )then
              call jsplt2( lgsigma )
            else
              call jsplt2( lrhohal )
            end if
          else
            if( disc( icmp ) )then
              call jsplt2( gmassd )
            else
              if( logm )then
                call jsplt2( lgmassh )
              else
                call jsplt2( gmassh )
              end if
            end if
          end if
          call jsdash( 0, 0, 0, 0 )
        end if
      end if
c switch to next run or original if done
      krun = krun + 1
      if( krun .le. nruns )then
        if( krun .eq. 2 )t = ts
        call cmprun( krun )
c allow for possible change of time step
        if( t .ne. ts )then
          jstep = nint( t * real( jstep ) / ts )
          jstep = max( jstep, 1 )
          iend = nint( t * real( iend ) / ts )
        end if
        call rewnd
        call nextrec( type, ifail )
        go to 5
      else
        call cmprun( 1 )
      end if
      return
      end

      real*8 function lgmassh( a )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
      real*8 a, r, gmassh, m
      r = 10.**a
      m = gmassh( r )
      lgmassh = log10( m )
      return
      end
