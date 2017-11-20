      subroutine spiral( is, ac, old )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the spiral field following Minchev & Famaey (2010).
c
c calling arguments
      integer is
      logical old
      real ac( 5 )
c      equivalence ( ac( 4 ), phi ), ( ac( 5 ), tau )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/model.f'
c
c local variables
      integer jstep
      real alpha, ar, arg, at, phi, r, scoeff, trot, x, y
      save jstep, scoeff
      include 'inc/pi.f'
c choose alpha=-2m, so pitch angle is 26 deg and negative to make it trail
      data alpha, jstep / -8., -1 /
c
      if( istep .ne. jstep )then
        if( .not. twod )call crash( 'SPIRAL', 'Not 2-D' )
        jstep = istep
c linearly growing disturbance
        trot = 2. * pi / omegsp
        arg = ts * real( istep ) / ( 3. * trot )
        arg = min( arg, 1. )
c        if( arg .lt. 1. )then
c          arg = arg**2 * ( 3. - 2. * arg )
c        else
c          arg = 1
c        end if
c scaled amplitude
        scoeff = arg * epspi
      end if
c
      if( old )then
        x = oldc( 1, is ) / lscale
        y = oldc( 2, is ) / lscale
      else
        x = newc( 1, is ) / lscale
        y = newc( 2, is ) / lscale
      end if
c evaluate expression
      r = sqrt( x**2 + y**2 )
      if( r .gt. 0. )then
        phi = atan2( y, x )
        arg = alpha * log( r ) - 4. * ( phi - bphase )
c acceleration components
        ar = scoeff * alpha * sin( arg ) / r
        at = -scoeff * 4. * sin( arg ) / r
        ac( 1 ) = ( ar * x - at * y ) / r
        ac( 2 ) = ( ar * y + at * x ) / r
c potential
        ac( 4 ) = -scoeff * sin( arg )
      else
        ac( 1 ) = 0
        ac( 2 ) = 0
        ac( 4 ) = 0
      end if
      ac( 3 ) = 0
c compute net torque produced by the spiral
      ac( 5 ) = x * ac( 2 ) - y * ac( 1 )
      return
      end
