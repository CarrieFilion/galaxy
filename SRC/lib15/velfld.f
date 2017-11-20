      subroutine velfld
c  Copyright (C) 201c, Jerry Sellwood
      use aarrays
      implicit none
c part of the analysis software
c   draws vectors of the mean stream speed of disk particles in a frame
c   that rotates with the bar
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/model.f'
c
c local variables
      integer ia, ifail, ir, j, ja, jstep, k, kstep, llstp, m
      real cth, dmax, omegap, phase, r, sth, th, vscale, vx, vy, x, y
      include 'inc/pi.f'
c
      call rewnd
      print *, 'Enter pattern speed'
      read *, omegap
    3 print *, 'Enter first step no, interval and last step'
      read *, jstep, kstep, llstp
c search for data
    1 call nextrec( 'LGSP', ifail )
      if( ifail .ne. 0 .or. istep .gt. llstp )return
      if( istep .lt. jstep )go to 1
c find bar major axis, given by phase of m = 2, p = 0 component
      m = 2
      j = ( mr / 2 ) * ma + m
      k = j + ma * mr
      phase = atan2( wres( k ), wres( j ) )
      if( phase .lt. 0. )phase = phase + 2. * pi
      phase = phase / real( m )
c now read in velocity field
      call nextrec( 'VLFD', ifail )
      if( istep .ne. jstep )print *, 'Relevant velocity field not found'
c linear scale
      dmax = real( mr ) / ( lscale * drfac )
      dmax = 3
c find bin closest to bar axis
      ja = phase * real( ma ) / ( 2. * pi )
      ja = ja + 1
c      print *, 'phase of bar  = ', phase * 180 / pi,
c     +        ' degrees,  or bin no', ja
c inclination of bar to nearest bin axis
c      pa = phase - 2. * pi * ( real( ja ) - .5 ) / real( ma )
c scale size of velocity vectors
      vscale = .5
c next frame
      call jspage
      call jssize( .15, 1., .1, .95 )
      call jscale( -dmax, dmax, -dmax, dmax )
c write heading
      call jsbldt( 'run no' )
      call jsbldi( irun, 5 )
      call jsbldt( 'time' )
      call jsbldf( time, 6, 2 )
      call jswrit( -.8 * dmax, 1.2 * dmax )
c draw and mark axes
      call jsaxis( 'x', ' ', 1 )
      call jsaxis( 'y', ' ', 1 )
c work over radii
      k = -5
      do 2 ir = 1, mr
      r = real( ir ) / ( lscale * drfac )
      if( r.ge.dmax )go to 3
c work over azimuthal bins
      do 5 ia = 1, ma
      k = k + 5
      if( wres( k + 1 ).lt.5. )go to 5
c angle to bar axis
      th = 2. * pi * ( real( ia ) - .5 ) / real( ma ) - phase
      cth = cos( th )
      sth = sin( th )
c cartesian components of streaming speed
      vx = .5 * vscale * ( cth * wres( k + 4 ) -
     +                     sth * ( wres( k + 2 ) - r * omegap ) )
      vy = .5 * vscale * ( sth * wres( k + 4 ) +
     +                     cth * ( wres( k + 2 ) - r * omegap ) )
      x = r * cth
      y = r * sth
      call jsarrw( x - vx, y - vy, x + vx, y + vy )
    5 continue
    2 continue
      jstep = jstep + llstp
      go to 3
      end
