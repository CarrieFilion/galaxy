      subroutine barquad( is, ac, old )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes an approximation of the bar quadrupole field following Weinberg
c   (1985) or Hernquist & Weinberg (1992).  Weinberg's (1985) prescription
c   for the parameter values is implemented in subroutine qudini_
c   returned values are in natural units
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
      logical homog
      real c2p, eta, gtime, r, scoeff, s2p, x, xi, y, z
      save c2p, jstep, scoeff, s2p
c      real*8 errfn, t
      parameter ( gtime = 10 )
      data jstep / -1 /, homog / .true. /
c
      if( istep .ne. jstep )then
        jstep = istep
c growing bar
        eta = ts * real( istep ) / gtime
        if( eta .lt. 1. )then
          eta = eta**2 * ( 3. - 2. * eta )
        else
          eta = 1
        end if
c growing bar - 42.886 of my time units = 1 in WK06
c        t = ts * real( istep ) / 42.886
c        t = 4. * t - 2.
c        i = 0
c        eta = .5 * ( 1. + errfn( t, i ) )
c scaled amplitude
        scoeff = eta * mptbr * alpha2 / abar**3
c trig factors
        c2p = cos( 2. * bphase )
        s2p = sin( 2. * bphase )
      end if
c
      if( old )then
        x = oldc( 1, is ) / lscale
        y = oldc( 2, is ) / lscale
        z = oldc( 3, is ) / lscale
      else
        x = newc( 1, is ) / lscale
        y = newc( 2, is ) / lscale
        z = newc( 3, is ) / lscale
      end if
      if( homog )then
c expression for a homogeneous ellipsoid given by Weinberg (1985)
        r = sqrt( x**2 + y**2 + z**2 )
        eta = r / ( beta2 * abar )
        xi = ( ( x**2 - y**2 ) * c2p + 2. * x * y * s2p ) / abar**2
        if( r .gt. 0. )then
c acceleration components
          ac( 1 ) = scoeff *
     +         ( 2. * ( 1. + eta**5 ) * ( x * c2p + y * s2p ) -
     +           5. * xi * x  * abar * eta**4 / ( r * beta2 ) ) /
     +                                                ( 1. + eta**5 )**2
          ac( 2 ) = scoeff *
     +         ( 2. * ( 1. + eta**5 ) * ( x * s2p - y * c2p ) -
     +           5. * xi * y  * abar * eta**4 / ( r * beta2 ) ) /
     +                                                ( 1. + eta**5 )**2
          ac( 3 ) = -scoeff * 5. * xi * z  * abar * eta**4 /
     +                                ( r * beta2 * ( 1. + eta**5 )**2 )
        else
          ac( 1 ) = 0
          ac( 2 ) = 0
          ac( 3 ) = 0
        end if
c potential
        ac( 4 ) = -scoeff * xi * abar**2 / ( 1. + eta**5 )
      else
c computes an approximation to the quadrupole field of an inhomogeneous
c   bar following Hernquist & Weinberg (1992)
        eta = beta2**2 + ( x**2 + y**2 + z**2 ) / abar**2
        xi = ( ( x**2 - y**2 ) * c2p + 2. * x * y * s2p ) / abar**2
c acceleration components
        ac( 1 ) = scoeff *
     +     ( 2. * eta * ( x * c2p + y * s2p ) - 5. * xi * x ) / eta**3.5
        ac( 2 ) = scoeff *
     +     ( 2. * eta * ( x * s2p - y * c2p ) - 5. * xi * y ) / eta**3.5
        ac( 3 ) = -scoeff * 5. * z * xi / eta**3.5
c potential
        ac( 4 ) = -scoeff * xi * abar**2 / eta**2.5
c old wrong expression
c      rc2 = x**2 + y**2
c      r2 = rc2 + z**2
c      d = betab**2 + rc2 / abar**2
c      e = ( x**2 - y**2 ) * c2p + 2. * x * y * s2p
cc acceleration components
c      ac( 1 ) = scoeff *
c     +          ( 2. * d * r2 * ( x * e + rc2 * ( x * c2p + y * s2p ) ) -
c     +          rc2 * e * x * ( 5. * r2 / abar**2 + 2. * d ) ) /
c     +                                                ( r2**2 * d**3.5 )
c      ac( 2 ) = scoeff *
c     +          ( 2. * d * r2 * ( y * e + rc2 * ( x * s2p - y * c2p ) ) -
c     +          rc2 * e * y * ( 5. * r2 / abar**2 + 2. * d ) ) /
c     +                                                ( r2**2 * d**3.5 )
c      ac( 3 ) = -scoeff * 2. * z * rc2 * e / ( r2**2 * d**2.5 )
cc potential
c      ac( 4 ) = -scoeff * rc2 * e / ( r2 * d**2.5 )
      end if
c compute net torque produced by the bar
      ac( 5 ) = x * ac( 2 ) - y * ac( 1 )
      return
      end
