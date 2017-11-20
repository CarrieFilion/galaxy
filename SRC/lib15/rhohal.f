      real*8 function rhohal( r )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c density of a spherical mass distribution
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      external rhalfn
c      real*8 aitken2
      real*8 deriv2, linint2, poveda, rhoisp, schmidt, sphrho, splint2
c
c local variables
      real*8 ai2, cc, haloc, hmass, hrad, r1
      include 'inc/pi.f'
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'RHOHAL',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        call crash( 'RHOHAL',
     +                     'Halo function called for a disc component' )
      else
        hmass = cmpmas( icmp )
        hrad = rscale( icmp )
      end if
c uniform sphere
      if( ctype( icmp ) .eq. 'UNIS' )then
        r1 = r / hrad
        rhohal = 0
        if( r1 .lt. 1.d0 )rhohal = .75 * hmass / ( pi * hrad**3 )
c Plummer sphere
      else if( ctype( icmp ) .eq. 'PLUM' )then
        r1 = r / hrad
        rhohal = hmass * .75 / ( pi * hrad**3 * ( 1. + r1**2 )**2.5 )
c cored isothermal
      else if( ctype( icmp ) .eq. 'ISOT' )then
        r1 = r / hrad
        rhohal = cmpmas( icmp ) *
     +        ( r1**2 + 3. ) / ( 4. * pi * hrad**3 * ( 1. + r1**2 )**2 )
c Jaffe model
      else if( ctype( icmp ) .eq. 'JAFF' )then
        rhohal = hmass * hrad / ( 4. * pi * r**2 * ( hrad + r )**2 )
c Kepler
      else if( ctype( icmp ) .eq. 'KEPL' )then
        rhohal = 0
c r**1/4 law
      else if( ctype( icmp ) .eq. 'RQUA' )then
        rhohal = hmass * poveda( r )
c modified Schmidt law
      else if( ctype( icmp ) .eq. 'SCHM' )then
        rhohal = hmass * schmidt( r )
c DFIT model
      else if( ctype( icmp ) .eq. 'DFIT' )then
       rhohal = sphrho( r, 0.d0, 0.d0 )
c King model or spherical polytrope - already scaled
      else if( ( ctype( icmp ) .eq. 'KING' ) .or.
     +         ( ctype( icmp ) .eq. 'POLY' ) )then
        rhohal = rhoisp( r )
c Hernquist model
      else if( ctype( icmp ) .eq. 'HERN' )then
        rhohal = hmass * hrad /
     +                ( 2. * pi * ( r + 1.d-10 ) * ( hrad + r )**3 )
c adiabatically (de-)compressed halo or VK03 model A_1
      else if( ( ctype( icmp ) .eq. 'ADIA' ) .or.
     +         ( ctype( icmp ) .eq. 'VKA1' ) )then
        if( r .gt. arad( nradad ) )then
          rhohal = 0
        else if( r0core )then
c simple interpolation, with extrapolation to zero radius
c          rhohal = aitken2( arad, rhon, nradad, r )
c spline interpolation - assume knots and coeffs are already set
          r1 = max( r, arad( 1 ) )
          rhohal = splint2( arad, rhon( 1, icomp ), nradad, r1,
     +                  rlamda( 1, icomp ), rhonc( 1, icomp ), .false. )
        else if( r1cusp )then
          if( r .lt. arad( 2 ) )then
c linear extrapolation to zero radius
            r1 = ( rhon( 2, icomp ) - rhon( 1, icomp ) ) /
     +                       ( log10( arad( 2 ) ) - log10( arad( 1 ) ) )
            rhohal = rhon( 1, icomp ) +
     +                          r1 * ( log10( r ) - log10( arad( 1 ) ) )
          else
c linear interpolation in the log
            rhohal = linint2( arad, rhon( 1, icomp ), nradad, r )
          end if
          rhohal = 10.d0**rhohal
        else
          call crash( 'RHOHAL', 'Unknown compressed halo form' )
        end if
c Henon's spherical isochrone
      else if( ctype( icmp ) .eq. 'ISOC' )then
        r1 = sqrt( 1. + r * r )
        rhohal = hmass * ( 1. + 2. * r1 ) /
     +                       ( 4. * pi * r1**3 * ( 1. + r1 )**2 )
c NFW profile
      else if( ctype( icmp ) .eq. 'NFW ' )then
        r1 = r / hrad + 1.d-15
        rhohal = hmass / ( 4. * pi * hrad**3 * r1 * ( 1. + r1 )**2 )
c singular isothermal sphere
      else if( ctype( icmp ) .eq. 'SISP' )then
        r1 = r + 1.d-10
        rhohal = hmass / ( 4. * pi * r1**2 )
c Hernquist's cored isothermal with a Gaussian cutoff
      else if( ctype( icmp ) .eq. 'ISOG' )then
        r1 = r / hrad
        cc = cmppar( 1, icmp )
        rhohal = hmass * cmppar( 2, icmp ) * exp( -r1**2 ) /
     +                  ( 2.d0 * pi**1.5 * hrad**3 * ( r1**2 + cc**2 ) )
c Einasto profile
      else if( ctype( icmp ) .eq. 'EINA' )then
        haloc = dfcns( 3, icmp )
        ai2 = 2.d0 / haloc
        r1 = r / hrad
        rhohal = hmass * exp( -ai2 * ( r1**haloc - 1. ) )
        rhohal = rhohal / ( 16.d0 * pi * hrad**3 )
c Burkert profile
      else if( ctype( icmp ) .eq. 'BURK' )then
        rhohal = hmass /
     +                 ( 4.d0 * pi * ( r + hrad ) * ( r**2 + hrad**2 ) )
c cubic halo - density tapers to zero at hrad
      else if( ctype( icmp ) .eq. 'CUBI' )then
        r1 = abs( r ) / hrad
        if( r1 .lt. 1.d0 )then
          rhohal = 15. * hmass * ( 2. * r1**3 - 3. * r1**2 + 1. ) /
     +                                             ( 4. * pi * hrad**3 )
        else
          rhohal = 0
        end if
      else
c general formula for spherical mass distribution
        r1 = max( r, 1.d-4 )
        rhohal = deriv2( r1, rhalfn ) / ( 4. * pi * r1**2 )
        if( rhohal .lt. 0.d0 )then
          print *, 'rhohal -ve for r =', sngl( r ), sngl( rhohal )
          call crash( 'RHOHAL', 'Halo density -ve' )
        end if
      end if
      return
      end
