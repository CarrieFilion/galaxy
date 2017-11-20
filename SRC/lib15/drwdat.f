      subroutine drwdat( mact )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to draw Fourier coefficients.  Used in data analysis software only.
c   Called from SELECT in MODEFIT or ANALYS
      use aarrays
      implicit none
c
c calling arguments
      integer mact
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
c external
      real roundup
c      real grofu
c
c local variables
      integer ip, it, n
      real x, xmax, xmin, y, ymax, ymin
c
c set y-axis scaling
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
c set x-axis scale
      xmin = 0
      xmax = tme( nt )
c write heading
      call jspage
      call jssize( 0., .75, 0., 1. )
      call jscale( 0., 1., 0., 1. )
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      call jsbldt( 'm =')
      call jsbldi( mact, 2 )
      call jsbldt( 'Time interval' )
      x = ( tme( kt ) - tme( jt ) ) / real( kt - jt )
      call jsbldf( x, 6, 1 )
      call jswrit( .1, .93 )
c set up frame
      call jssize( .1, .9, .1, .85 )
      call jscale( xmin, xmax, ymin, ymax )
c draw and mark axes
      call jsaxis( 'X', 'Time', 1 )
      call jsaxis( 'Y', 'ln(A)', 1 )
c plot amplitude as a function of time
      do ip = jp, kp
        n = ip - jp + 1
        do it = jt, kt
          y = sqrt( sdata( 1, n )**2 + sdata( 2, n )**2 )
          y = max( y, .1e-10 )
          y = log( y )
          y = min( y, ymax )
          x = tme( it )
          if( it .eq. jt )call jsmove( x, y )
          call jsline( x, y )
          n = n + kp - jp + 1
        end do
      end do
      return
      end
