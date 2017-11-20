      real*8 function Phihal( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns gravitational potential from a spherical mass component at the
c   requested radius
c   may use QUADPAK routine DQAGI
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
      include 'inc/model.f'
c
c externals
      real*8 frhalo, Gam1o1, gammaP, gmassh, isopsi
      real*8 Phidsc, Phisph, quad_inf, splint2
c
c local arrays
      integer ntable, iuseph
      real*8 maxr
      parameter ( ntable = 501, maxr = 20. )
      real*8, allocatable :: radtab(:), Phitab(:), ctab(:), lamtab(:)
      save iuseph, ctab, lamtab, Phitab, radtab
c
c scratch arrays
      real*8, allocatable :: absc( : ), wt( : )
c
c local variables
      integer i, ier, jcmp, npts
      real*8 ai, epsa, epsr, haloc, hmass, hrad, r1, r2, xg
      parameter ( npts = 32 )
      data iuseph / 0 /
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'PHIHAL',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        call crash( 'PHIHAL', 'Called for a disc component' )
      else
        hmass = cmpmas( icmp )
        hrad = rscale( icmp )
        haloc = dfcns( 3, icmp )
      end if
c
      r1 = r / hrad
c fixed rotation curve models - assume halo cuts off sharply at rmax
      if( fixed )then
        Phihal = rmax * frhalo( rmax )
        if( r .ge. rmax )then
          Phihal = Phihal * rmax / r
        else
c add integral of radial force from r to rmax
          r2 = r
          allocate ( absc( npts ) )
          allocate ( wt( npts ) )
          call GLqtab( r2, rmax, npts, wt, absc )
          do i = 1, npts
            r2 = absc( i )
            Phihal = Phihal + wt( i ) * frhalo( r2 )
          end do
        end if
c frozen disc
      else if( ctype( icmp ) .eq. 'ASDI' )then
        jcmp = icmp
        do icmp = 1, ncmp
          if( disc( icmp ) .and. ( cmpmas( icmp ) .gt. 0. ) )Phihal =
     +                     cmpmas( jcmp ) * Phidsc( r ) / cmpmas( icmp )
        end do
        icmp = jcmp
c uniform sphere
      else if( ctype( icmp ) .eq. 'UNIS' )then
        if( r1 .lt. 1. )then
          Phihal = hmass * ( .5 * r1 * r1 - 1.5 ) / hrad
        else
          Phihal = -hmass / r
        end if
c Plummer sphere
      else if( ctype( icmp ) .eq. 'PLUM' )then
        Phihal = -hmass / sqrt( hrad * hrad + r * r )
c cored isothermal halo
      else if( ctype( icmp ) .eq. 'ISOT' )then
        r1 = r / hrad
        Phihal = .5 * hmass * log( 1. + r1**2 ) / hrad
c Jaffe model
      else if( ctype( icmp ) .eq. 'JAFF' )then
        r1 = max( r, 1.d-12 )
        Phihal = hmass * log( r1 / ( hrad + r ) ) / hrad
c Kepler
      else if( ctype( icmp ) .eq. 'KEPL' )then
        r1 = max( r, 1.d-12 )
        Phihal = -hmass / r1
c r**1/4 law
      else if( ctype( icmp ) .eq. 'RQUA' )then
        epsa = 1.d-8
        epsr = 1.d-8
        Phihal = quad_inf( frhalo, r, 1, epsa, epsr, ier )
        if( ier .gt. 0 )then
          print *, 'ier =', ier, ' from QUAD_INF'
          call crash( 'PHIHAL', 'QUADPACK error' )
        end if
c DFIT model - use Plummer sphere expression to provide a non-zero value
      else if( ctype( icmp ) .eq. 'DFIT' )then
        Phihal = -hmass / sqrt( hrad * hrad + r * r )
c TEDS model
      else if( ctype( icmp ) .eq. 'TEDS' )then
        r1 = hrad * haloc
        Phihal = .5 / sqrt( hrad * hrad + r * r ) +
     +          .09 / sqrt( r1 * r1 + r * r )
        Phihal = -hmass * Phihal / .59
c King model
      else if( ctype( icmp ) .eq. 'KING' )then
        Phihal = -hmass / rtidal - isopsi( r )
c polytropic sphere
      else if( ctype( icmp ) .eq. 'POLY' )then
        Phihal = -hmass / rscale( icmp ) - isopsi( r )
c Kuz'min-Kutuzov oblate sphere - potential in midplane
      else if( ctype( icmp ) .eq. 'KKSP' )then
        r1 = 0
        Phihal = hmass * Phisph( r, r1 )
c Hernquist model
      else if( ctype( icmp ) .eq. 'HERN' )then
        Phihal = -hmass / ( hrad + r )
c UNKN model - use Plummer sphere expression to provide a non-zero value
      else if( ctype( icmp ) .eq. 'UNKN' )then
        Phihal = -hmass / sqrt( hrad * hrad + r * r )
c generalized Aguilar-Merritt model
      else if( ctype( icmp ) .eq. 'AGME' )then
        if( r1 .lt. 1. )then
          Phihal = -hmass *
     +      ( 1. + ( 1. - r1**( 2. + haloc ) ) / ( 2. + haloc ) ) / hrad
        else
          Phihal = -hmass / r
        end if
c adiabatically (de-)compressed halo - potential of separate compments unknown
      else if( ctype( icmp ) .eq. 'ADIA' )then
        call crash( 'PHIHAL', 'Logic error' )
c Henon's spherical isochrone
      else if( ctype( icmp ) .eq. 'ISOC' )then
        Phihal = -hmass / ( hrad + sqrt( hrad * hrad + r * r ) )
c NFW profile
      else if( ctype( icmp ) .eq. 'NFW ' )then
        r1 = r / hrad
        if( abs( r1 ) .gt. 1.d-5 )then
          Phihal = -hmass * log( 1.d0 + r1 ) / r
        else
c expansion for small arguments
          Phihal = -hmass *
     + ( 1.d0  - 0.5 * r + r**2 / 3. - 0.25 * r**3 + 0.2 * r**4 ) / hrad
        end if
c singular isothermal sphere
      else if( ctype( icmp ) .eq. 'SISP' )then
        r1 = 1.e-20
        r1 = max( r, r1 )
        Phihal = hmass * log( r1 )
c Miyamoto-Nagai spheroid - potential in midplane
      else if( ctype( icmp ) .eq. 'MYNG' )then
        r1 = 0
        Phihal = hmass * Phisph( r, r1 )
c Einasto halo - from Cardone et al (2005, MNRAS, 358, 1325)
      else if( ctype( icmp ) .eq. 'EINA' )then
        r1 = r / hrad
        ai = 1.d0 / haloc
        xg = 2.d0 * ai * r1**haloc
        ier = 0
        phihal = ( 1.d0 - gammaP( 2.d0 * ai, xg, ier ) ) *
     +                                          ( 2.d0 * ai )**ai / hrad
        phihal = phihal * Gam1o1( 2.d0 * ai, 3.d0 * ai )
        phihal = -hmass * mnorm( icmp ) * phihal -
     +            gmassh( r ) / ( r + 1.d-20 )
c cubic halo - density tapers to zero at hrad
      else if( ctype( icmp ) .eq. 'CUBI' )then
        r1 = abs( r ) / hrad
        if( r1 .lt. 1.d0 )then
          ai = 9. - r1**2 * ( 4. * r1**3 - 9. * r1**2 + 10. )
          phihal = -hmass * ai / ( 4. * hrad )
        else
          phihal = -hmass / r
        end if
c integrate radial force
      else if( .not. sphrod( icmp ) )then
        if( iuseph .eq. 0 )then
          if( master )print *, 'PHIHAL: Building a table'
          allocate ( radtab( ntable ) )
          allocate ( Phitab( ntable ) )
c numerical integral of halo radial force
          do i = 1, ntable
            r1 = maxr * ( real( i - 1 ) / real( ntable - 1 ) )**2
            epsa = 1.d-10
            epsr = 1.d-10
            Phihal = quad_inf( frhalo, r1, 1, epsa, epsr, ier )
            if( ier .gt. 0 )then
              print *, 'ier =', ier, ' from QUAD_INF'
              call crash( 'PHIHAL', 'QUADPACK error' )
            end if
            radtab( i ) = r1
            Phitab( i ) = Phihal
          end do
          allocate ( ctab( ntable + 4 ) )
          allocate ( lamtab( ntable + 4 ) )
          r1 = .5 * r1
          Phihal = splint2( radtab, Phitab, ntable, r1,
     +                      lamtab, ctab, .true. )
          if( master )print *, 'PHIHAL: Table ready'
          iuseph = 1
        end if
        if( r .le. maxr )then
          Phihal = splint2( radtab, Phitab, ntable, r,
     +                      lamtab, ctab, .false. )
        else
          epsa = 1.d-8
          epsr = 1.d-8
          Phihal = quad_inf( frhalo, r, 1, epsa, epsr, ier )
          if( ier .gt. 0 )then
            print *, 'ier =', ier, ' from QUAD_INF'
            call crash( 'PHIHAL', 'QUADPACK error' )
          end if
        end if
      else
c none of the above
        call crash( 'PHIHAL', 'Unknown halo type' )
      end if
      return
      end
