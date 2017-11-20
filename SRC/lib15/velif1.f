      real*8 function velif1( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the integral of a, possibly velocity weighted, DF over the allowed
c   range of azimuthal velocities at fixed radius and radial velocity.  It
c   works for all disk and spheroidal models with two-integral DFs
c   uses QUADPACK routine DQNG
c
c calling argument
      real*8 x
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
c externals
      external velif2
      real*8 quad_nad
c
c local variables
      integer ier
      real*8 epsa, epsr, vmax, vmin
c
c Gauss-Legendre quadrature over allowed range of v
      u = x
c maximum tangential velocity - if circular orbits at edge are included
      if( emaxe .gt. phimax )then
        vmax = rmax * sqrt( ( 2. * ( phimax - pot ) - u * u )
     +                                 / ( rmax * rmax - rad * rad ) )
c additional check (in case of small denominator above)
        if( ( rad .gt. 0.d0 ) .and. ( emaxe .gt. phimax ) )vmax =
     +   min( vmax, ( rmax / rad ) * sqrt( 2. * ( emaxe - phimax ) ) )
      else
c if just a simple energy cut-off
        vmax = sqrt( 2. * ( phimax - pot ) - u * u )
      end if
c minumum value
      vmin = -vmax
      if( .not. retract )then
        if( rad .eq. 0.d0 )then
          if( Lzrmn( icmp ) .eq. 0. )vmin = 0.
        else
          vmin = max( -vmax, Lzrmn( icmp ) / rad )
        end if
      end if
c set accuracy criteria impossibly high
      epsa = 1.e-8
      epsr = epsa
      velif1 = quad_nad( velif2, vmin, vmax, epsa, epsr, ier )
c ier = 1 means requested accuracy was not achieved
      if( ier .gt. 1 )then
        print *, 'ier =', ier, ' from QUAD_NAD'
        call crash( 'VELIF1', 'QUADPACK error' )
      end if
c double result if retrograde stars were omitted for spheroidal functions
      if( ( .not. disc( icmp ) ) .and.
     +    ( vmin .eq. 0. ) )velif1 = 2 * velif1
      return
      end
