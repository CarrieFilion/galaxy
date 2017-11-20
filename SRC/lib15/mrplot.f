      subroutine mrplot
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot the time variation of the radii containing various fractions
c   of the particles
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
      include 'inc/comprns.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      logical gtlogl
      real grofu, roundup
c
c local arrays
      integer mms
      parameter ( mms = 20 )
      real fm( mms ), xk( 2 ), yk( 2 )
      real, allocatable :: rad( :, : )
      real, allocatable :: tm( : )
c
c local variables
      character*4 type
      integer i, ifail, ip, isize, istp, j, krun, l, nms, nstride
      logical logy
      real a, area, b, r1, r2, xmax, ymax, ymin
c
c select type of data
      call datype( type )
c determine which population to plot
      call selpop( ip )
      icmp = ip
c get first record
      call rewnd
      call nextrec( type, ifail )
      if( ifail .ne. 0 )then
        print *, 'No data of the requested type found'
        return
      end if
c choose type of scale - log scale not allowed for danl data
      if( danl )then
        logy = .false.
      else
        logy = gtlogl( 'Logarithmic scale in radius' )
      end if
c select mass fractions
      a = fmass( ip ) * cmpmas( ip )
      if( danl )then
        nms = mms
      else
        a = nsp( icmp )
        if( logy )then
          b = a / real( ma )
          b = log10( b ) / log10( 2. )
          nms = nint( b )
          nms = min( nms, mms )
          nstride = ma
        else
          nms = mms
        end if
      end if
      if( logy )then
        do i = 1, nms
          fm( i ) = ma * 2**( i - 1 )
        end do
      else
        do i = 1, nms
          fm( i ) = a * ( real( i ) - .5 ) / real( nms )
        end do
      end if
c label plot
      call jspage
      call jssize( 0., 1., 0., 1. )
      call jscale( 0., 1., 0., 1. )
      if( nruns .eq. 1 )then
        call jsbldt( 'Run no' )
        call jsbldi( irun, 5 )
      end if
      call jsbldt( ' Mass fractions in component' )
      call jsbldi( icmp, 1 )
      call jswrit( .2, .95 )
c count number of records and allocate space
      istp = 0
      do while ( ifail .eq. 0 )
        istp = istp + 1
        call nextrec( type, ifail )
      end do
      isize = istp
      allocate ( rad( isize, mms ) )
      allocate ( tm( isize ) )
c
      ymin = 1000000
      ymax = 0
c work through file
      krun = 1
    1 istp = 0
      call rewnd
      call nextrec( type, ifail )
      do while ( ifail .eq. 0 )
        istp = istp + 1
        tm( istp ) = time
        do i = 1, nms
          rad( istp, i ) = 0
        end do
        if( danl )then
c convert surface densities back to masses
          r1 = grofu( 0. ) / lscale
          a = 0
          j = 1
          do i = 1, mr
            r2 = i
            r2 = grofu( r2 - .5 ) / lscale
c .5 * alpha = pi / na
            area = .5 * alpha * ( r2 * r2 - r1 * r1 )
            l = ( i - 1 ) * ma + 1
            wres( l ) = wres( l ) * area
c locate desired radii
            if( ( a + wres( l ) .gt. fm( j ) ) .and.
     +          ( j .le. nms ) )then
              rad( istp, j ) = r1 +
     +                         ( fm( j ) - a ) * ( r2 - r1 ) / wres( l )
              j = j + 1
            end if
            a = a + wres( l )
            r1 = r2
          end do
          if( krun .eq. 1 )then
            ymin = min( ymin, rad( istp, 1 ) )
            ymax = max( ymax, rad( istp, nms ) )
          end if
        else if( rhor .or. sigr )then
c rhor data
          r2 = 0
          b = 0
          j = 1
          do i = 1, mr, 2
            r1 = r2
            r2 = 2 * wres( i ) - r1
            a = b
            b = .5 * real( ma * ( i + 1 ) )
c locate desired radii
            if( ( b .gt. fm( j ) ) .and. ( j .le. nms ) )then
              rad( istp, j ) = r1 +
     +                         ( fm( j ) - a ) * ( r2 - r1 ) / ( b - a )
              j = j + 1
            end if
          end do
          if( krun .eq. 1 )then
            ymin = min( ymin, rad( istp, 1 ) )
            ymax = max( ymax, rad( istp, nms ) )
          end if
        end if
c get next record
        call nextrec( type, ifail )
        if( istp .eq. isize )ifail = 1
      end do
c convert to log scale
      if( logy )then
        do i = 1, istp
          do j = 1, nms
            if( rad( i, j ) .gt. 0. )then
              rad( i, j ) = log10( rad( i, j ) )
            else
              rad( i, j ) = 10
            end if
          end do
        end do
        if( krun .eq. 1 )then
          ymax = max( -1., log10( ymax ) )
          ymin = min( -2.5, log10( ymin ) )
        end if
      end if
c set up plot and draw axes
      if( krun .eq. 1 )then
        xmax = roundup( .99 * tm( istp ) )
        if( logy )then
          if( ymin .gt. 0. )then
            ymin = roundup( .5 * ymin )
          else
            ymin = -roundup( -1.01 * ymin )
          end if
          if( ymax .gt. 0 )ymax = roundup( .99 * ymax )
        else
          ymax = roundup( .99 * ymax )
c          ymax = max( ymax, rgrid( ngrid ) / lscale )
          ymin = 0
        end if
        call jssize( .15, .95, .15, .85 )
        call jscale( 0., xmax, ymin, ymax )
        call jsaxis( 'X', 'Time', 4 )
        call jsbldt( ctype( icmp ) )
        if( disc( icmp ) )then
          call jsbldt( 'disk:' )
        else
          call jsbldt( 'bulge/halo:' )
        end if
        if( rhor )call jsbldt( 'RHOR' )
        if( sigr )call jsbldt( 'SIGR' )
        if( danl )call jsbldt( 'DANL' )
        call jsbldt( 'data' )
        call jswrit( .2 * xmax, 1.05 * ymax - 0.05 * ymin )
        if( logy )then
          call jsaxis( 'Y', 'log_{10}r', 1 )
          if( rhor )then
            call jsbldt( '1st mass fraction = 1 /' )
            i = nsp( icmp ) / nstride
            a = log10( real( i ) )
            j = a + 1.01
            call jsbldi( i, j )
            call jsbldt( 'or' )
            a = log10( real( nstride ) )
            j = a + 1.01
            call jsbldi( nstride, j )
            call jsbldt( 'particles' )
            call jschsz( .2 )
            call jswrit( .5 * xmax, .05 * ymax + .95 * ymin )
            call jschsz( .3 )
          end if
        else
          call jsaxis( 'Y', 'r', 1 )
        end if
      end if
c identify run
      if( nruns .gt. 1 )then
        xk( 1 ) = .85 * xmax
        xk( 2 ) = .88 * xmax
        yk( 1 ) = ( .95 - .05 * real( krun ) )
        yk( 1 ) = yk( 1 ) * ymax + ( 1. - yk( 1 ) ) * ymin
        yk( 2 ) = yk( 1 )
        call jscoin( xk, yk, 2, krun )
        call jsbldi( irun, 4 )
        yk( 1 ) = ( .95 - .05 * real( krun ) ) - .01
        yk( 1 ) = yk( 1 ) * ymax + ( 1. - yk( 1 ) ) * ymin
        call jswrit( xk( 2 ) + .01 * xmax, yk( 1 ) )
      end if
c plot data
      do i = 1, nms
        if( logy )then
          j = 0
          do while ( ( j .lt. istp ) .and. ( rad( j + 1, i ) .lt. 5. ) )
            j = j + 1
          end do
          j = min( j, istp )
        else
          j = 0
          do while ( ( j .lt. istp ) .and. ( rad( j + 1, i ) .lt. 5. ) )
            j = j + 1
          end do
          j = min( j, istp )
        end if
        call jscoin( tm, rad( 1, i ), j, krun )
      end do
c switch to next run or original if done
      krun = krun + 1
      if( krun .le. nruns )then
        call cmprun( krun )
        go to 1
      else
        call cmprun( 1 )
        krun = 1
      end if
      return
      end
