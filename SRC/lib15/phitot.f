      real*8 function Phitot( r )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c returns potential in the mid-plane due to both disc and halo mass compnts
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      real*8 frtot, Phidsc, Phiext, Phihal, Phispl, splint2
c, Phitab
c
c local arrays
      integer npts
      parameter ( npts = 32 )
      real*8 absc( npts ), wt( npts )
c
c local variables
      integer i, jcmp, k
      real*8 gm, haloc, hrad, r1, r2
c
c preserve icmp
      jcmp = icmp
      if( numcal )then
        if( lgfld )then
c numerical potential, value in mid-plane only
          r2 = min( r, arad( nradad ) )
          Phitot =
     +         splint2( arad, Phin, nradad, r2, plamda, phinc, .false. )
c add integral of radial force from r2 to r
          if( r2 .lt. r )then
            k = 4
            call GLqtab( r2, r, k, wt, absc )
            do i = 1, k
              r2 = absc( i )
              Phitot = Phitot - wt( i ) * frtot( r2 )
            end do
          end if
        else if( ldfit )then
c numerical potential, value in mid-plane only for spheroidal systems
          do icmp = 1, ncmp
            if( sphrod( icmp ) )Phitot = Phispl( r, 0.d0 )
          end do
        else
          call crash( 'PHITOT', 'Unrecognized numcal option' )
        end if
      else if( fixed )then
c fixed rotation curve model
        icmp = 0
        do i = 1, ncmp
          if( .not. disc( i ) )then
            haloc = dfcns( 3, i )
            hrad = rscale( i )
            icmp = i
          end if
        end do
        if( icmp .eq. 0 )call crash( 'PHITOT',
     +            'No halo cmp found for a fixed rotation curve model' )
c Fall-Efstathiou model
        if( ctype( icmp ) .eq. 'FE  ' )then
          Phitot = .5 * haloc**2 * log( ( r**2 + hrad**2 )
     +                                         / ( rmax**2 + hrad**2 ) )
c power law rotation curve - choose zero at r = 1
        else if( ctype( icmp ) .eq. 'POWE' )then
          r1 = r
          if( haloc .lt. .01 )r1 = max( r, 1.d-12 )
          if( abs( haloc ) .lt. .01 )then
            Phitot = log( r1 )
          else
            Phitot = ( r1**( 2. * haloc ) - 1. ) / ( 2. * haloc )
          end if
        else
c assume halo cuts off sharply at rmax
          Phitot = rmax * frtot( rmax )
          if( r .ge. rmax )then
            Phitot = Phitot * rmax / r
          else
c add integral of radial force from r to rmax
            r2 = r
            call GLqtab( r2, rmax, npts, wt, absc )
            do i = 1, npts
              r2 = absc( i )
              Phitot = Phitot + wt( i ) * frtot( r2 )
            end do
          end if
        end if
c case needed only at start of compress
      else if( snglph )then
        Phitot = Phihal( r )
c adiabatically compressed halo - common potential for all components
      else if( cmprssd )then
        if( r .gt. arad( nradad ) )then
          Phitot = Phin( nradad )
c rigid components may have non-zero density at these radii
          do irigid = 1, nrigid
            icmp = ircmp( irigid )
            Phitot = Phitot - Phiext( arad( nradad ) ) + Phiext( r )
          end do
c compressible components may have zero density at these radii
          gm = 0
          do i = 1, ncomp
            gm = gm + gmn( nradad, i )
          end do
          Phitot = Phitot + gm * ( 1. / arad( nradad ) - 1. / r )
        else if( r .lt. arad( 1 ) )then
c linear extrapolation to zero radius
          r1 = ( Phin( 2 ) - Phin( 1 ) ) / ( arad( 2 ) - arad( 1 ) )
          Phitot = Phin( 1 ) + r1 * ( r - arad( 1 ) )
        else
c spline interpolation - assume knots and coeffs are already set
          Phitot = splint2( arad, Phin, nradad, r, plamda, Phinc,
     +                      .false. )
        end if
      else
        Phitot = 0
c regular model
        do icmp = 1, ncmp
          if( disc( icmp ) )then
c disc component
            Phitot = Phitot + Phidsc( r )
          else
            if( ctype( icmp ) .eq. 'POWC' )then
c power law with a core
              call crash( 'PHITOT',
     +                        'Power-law option not programmed' )
c              if( abs( dfcns( 3, icmp ) ) .lt. .01 )then
c                Phitot = .5 * ( log( 1. + ( r / hrad )**2 ) - 5. ) / hrad
c              else
c                Phitot = ( 1. + ( r / hrad )**2 )**( -.5 * dfcns( 3, icmp ) )
c                if( dfcns( 3, icmp ) .lt. 0.d0 )
c     +          Phitot = ( rmax / hrad )**( -dfcns( 3, icmp ) ) - Phitot
c                Phitot = -Phitot / hrad
c              end if
            else
c halo component
              Phitot = Phitot + Phihal( r )
            end if
          end if
        end do
      end if
c restore icmp
      icmp = jcmp
      return
      end
