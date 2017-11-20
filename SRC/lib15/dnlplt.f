      subroutine dnlplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot non-axisymmetric density amplitudes at a fixed radius
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
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      real roundup, uofgr
      real*8 akcrit
c
c local arrays
      integer mstp
      parameter ( mstp = 1001 )
      real amp( mstp ), phase( mstp ), tm( mstp )
c
c local variables
      integer ifail, ir, istp, j, k, m, mact
      real a, r, x, xmax, ymax, ymin
      real*8 r2
c
      icmp = 1
c read in Fourier harmonic
    1 call choosm( 'DANL', mact, m )
      if( mact .le. 0 )return
      call jschsz( .25 )
      ir = 0
      do while ( ( ir .gt. nr( jgrid ) ) .or. ( ir .lt. 1 ) )
        call gtreal( 'Select radius or zero to skip', r )
        if( r .le. 0. )go to 1
        ir = 1. + uofgr( r * lscale )
      end do
c read down file
      call rewnd
      istp = 0
      call nextrec( 'DANL', ifail )
      do while ( ifail .eq. 0 )
        istp = istp + 1
        if( istp .gt. mstp )call crash( 'DNLPLT', 'Arrays too small' )
        j = ( ir - 1 ) * ma + m + 1
        k = j + mr * ma
        a = sqrt( wres( j ) * wres( j ) + wres( k ) * wres( k ) )
        if( a .gt. 0. )then
          amp( istp ) = log( a )
          phase( istp ) = atan2( wres( k ), wres( j ) ) / real( m )
        else
          amp( istp ) = -10
          phase( istp ) = 0
        end if
        tm( istp ) = time
        call nextrec( 'DANL', ifail )
      end do
c plot growth of selected harmonic
      call jspage
      call jssize( .1, .7, .1, .9 )
      call jscale( 0., 1., 0., 1. )
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      call jsbldt( ' m =' )
      call jsbldi( mact, 2 )
      call jsbldt( ' r =' )
      call jsbldf( r, 6, 2 )
c compute and output Toomre's X from original disk
      r2 = r
      x = r * akcrit( r2 ) / real( mact )
      call jsbldt( ' x =' )
      call jsbldf( x, 6, 2 )
      call jswrit( .1, .8 )
c reset window
      xmax = roundup( tm( istp ) )
      call jssize( .15, .95, 0.1, 0.7 )
      ymax = 0.
      ymin = -5.
      call jscale( 0., xmax, ymin, ymax )
c draw and label axes
      call jsaxis( 'y', 'ln(a)', 1 )
      call jsaxis( 'x', 'time', 1 )
c plot amplitudes
      call jsymbl( tm, amp, istp, 4 )
      go to 1
      end
