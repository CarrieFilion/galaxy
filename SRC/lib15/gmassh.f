      real*8 function gmassh( r )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c returns halo mass interior to r
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
      integer ntable, iusegm
      real*8 maxr
      parameter ( ntable = 1001, maxr = 6 )
      real*8, allocatable :: radtab(:), gmtab(:)
      save iusegm, gmtab, radtab
c
c externals
      real*8 aitken2, algrng2, bik0, bik1, gammaP, gmfunh
      real*8 gmisph, gmtabh, splint2, vcirc, vdisc
c
c local arrays
      real*8, allocatable :: absc( : ), wt( : )
c
c local variables
      integer i, ifail, j, jcmp, npts
      real*8 ag, haloc, hmass, hrad, r1, rm, rs, xg
      include 'inc/pi.f'
      data iusegm / 0 /
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'GMASSH',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        call crash( 'GMASSH', 'Called for a disc component' )
      else
        hmass = cmpmas( icmp )
        hrad = rscale( icmp )
        haloc = dfcns( 3, icmp )
      end if
c
      gmassh = -1
      rs = r / hrad
      if( r .le. 0. )then
        gmassh = 0
c uniform sphere
      else if( ctype( icmp ) .eq. 'UNIS' )then
        rs = min( rs, 1.d0 )
        gmassh = hmass * rs**3
c Plummer sphere - use this as a sph approx to a Miyamoto-Nagai spheroid
      else if( ( ctype( icmp ) .eq. 'PLUM' ) .or.
     +         ( ctype( icmp ) .eq. 'MYNG' ) )then
        gmassh = hmass * rs**3 / ( 1. + rs * rs )**1.5
c simple 'isothermal' halo
      else if( ctype( icmp ) .eq. 'ISOT' )then
        gmassh = hmass * rs**3 / ( 1. + rs * rs )
c Jaffe model
      else if( ctype( icmp ) .eq. 'JAFF' )then
        gmassh = hmass * rs / ( 1. + rs )
c r**1/4 or Burkert models - create a look-up table at first call
      else if( ctype( icmp ) .eq. 'RQUA' .or.
     +         ctype( icmp ) .eq. 'BURK' )then
        if( iusegm .eq. 0 )then
          if( master )print *,'GMASSH: Building a table'
          allocate ( radtab( ntable ) )
          allocate ( gmtab( ntable ) )
          npts = 32
          allocate ( absc( npts ) )
          allocate ( wt( npts ) )
          rm = rtrunc( icmp )
          rm = max( rm, maxr )
c non-adaptive Gauss-Legendre quadrature of halo density expression
          do i = 1, ntable
            rs = rm * real( i - 1 ) / real( ntable - 1 )
            r1 = 0
            call GLqtab( r1, rs, npts, wt, absc )
            gmassh = 0
            do j = 1, npts
              r1 = absc( j )
              gmassh = gmassh + wt( j ) * gmfunh( r1 )
            end do
            radtab( i ) = rs
            gmtab( i ) = 4. * pi *gmassh
          end do
          iusegm = 1
          rs = r / hrad
          if( master )print *,'GMASSH: Table ready'
        end if
        if( rs .le. maxr )then
          if( rs .gt. radtab( 2 ) )then
            gmassh = algrng2( radtab, gmtab, ntable, rs )
          else
            if( ctype( icmp ) .eq. 'BURK' )then
c to 2nd order, rho propto (1-rs), so Taylor series gives
              xg = radtab( 2 )
              gmassh = gmtab( 2 ) * ( rs / xg )**3 * ( 1. - .75 * rs ) /
     +                                               ( 1. - .75 * xg )
            else if( ctype( icmp ) .eq. 'RQUA' )then
              gmassh = gmtab( 2 ) * ( rs / radtab( 2 ) )**2
            else
              call crash( 'GMASSH', 'Need expression for small r' )
            end if
          end if
        else
          r1 = radtab( ntable )
          npts = 8
          allocate ( absc( npts ) )
          allocate ( wt( npts ) )
          call GLqtab( r1, r, npts, wt, absc )
          gmassh = 0
          do j = 1, npts
            r1 = absc( j )
            gmassh = gmassh + wt( j ) * gmfunh( r1 )
          end do
          gmassh = 4. * pi * gmassh + gmtab( ntable )
        end if
c TEDS model
      else if( ctype( icmp ) .eq. 'TEDS' )then
        gmassh = .5 * rs**3 / ( 1. + rs * rs )**1.5
        rs = rs * haloc
        gmassh = gmassh + .09 * rs**3 / ( 1. + rs * rs )**1.5
        gmassh = hmass * gmassh / .59
c King model or spherical polytrope
      else if( ( ctype( icmp ) .eq. 'KING' ) .or.
     +         ( ctype( icmp ) .eq. 'POLY' ) )then
        gmassh = gmisph( r )
c Hernquist model
      else if( ctype( icmp ) .eq. 'HERN' )then
        gmassh = hmass * rs**2 / ( 1. + rs )**2
c generalized Aguilar-Merritt model
      else if( ctype( icmp ) .eq. 'AGME' )then
        gmassh = hmass
        if( rs .lt. 1. )gmassh = hmass * rs**( 3. + haloc )
      else if( ( ctype( icmp ) .eq. 'ADIA' ) .or.
     +         ( ctype( icmp ) .eq. 'VKA1' ) )then
c adiabatically (de-)compressed halo or VK03 model A_1
        if( r .gt. arad( nradad ) )then
          gmassh = gmn( nradad, icomp )
        else if( r0core )then
          if( r .gt. arad( 2 ) )then
            gmassh = aitken2( arad, gmn( 1, icomp ), nradad, r )
          else
c assume a linear density gradient over 0 < r < rad( 2 )
            gmassh = pi * ( 4. * r**3 * rhon( 1, icomp ) / 3. +
     +                r**4 * ( rhon( 2, icomp ) - rhon( 1, icomp ) ) /
     +                                       ( arad( 2 ) - arad( 1 ) ) )
          end if
        else if( r1cusp )then
          if( r .gt. arad( 1 ) )then
c spline interpolation - assume knots and coeffs are already set
            gmassh = splint2( arad, gmn( 1, icomp ), nradad, r,
     +                        glamda( 1, icomp ), gmnc( 1, icomp ),
     +                        .false. )
          else
c assume an r^{-1} density gradient over 0 < r < rad( 1 )
            gmassh = 2. * pi * r**2 * rhon( 1, icomp ) * arad( 1 )
          end if
        else
          call crash( 'GMASSH', 'Unknown compressed halo form' )
        end if
c spherical isochrone
      else if( ctype( icmp ) .eq. 'ISOC' )then
        r1 = sqrt( hrad * hrad + r * r )
        gmassh = hmass * r**3 / ( r1 * ( hrad + r1 )**2 )
c NFW profile
      else if( ctype( icmp ) .eq. 'NFW' )then
        r1 = r / hrad
        if( r1 .gt. 1.d-5 )then
          gmassh = log( 1. + r1 ) - r1 / ( 1. + r1 )
        else
          gmassh = .5 * r1**2 - 2. * r1**3 / 3. + .75 * r1**4
        end if
        gmassh = hmass * gmassh
c singular isothermal sphere
      else if( ctype( icmp ) .eq. 'SISP' )then
        gmassh = hmass * r
c Einasto finite central density model
      else if( ctype( icmp ) .eq. 'EINA' )then
        r1 = r / hrad
        ag = 3.d0 / haloc
        xg = 2.d0 * r1**haloc / haloc
c incomplete gamma function = P(a;x) * Gamma(a)
        ifail = 0
        gmassh = hmass * mnorm( icmp ) * gammaP( ag, xg, ifail )
      else if( ctype( icmp ) .eq. 'MTAB' )then
c tabulated mass profiles
        gmassh = hmass * gmtabh( r )
c halo adopted by Donner-Thomasson - uses potential of a double exp disk!
      else if( ctype( icmp ) .eq. 'EXPH' )then
        r1 = .5 * abs( r / hrad )
        gmassh = 0.2 * r1**2 * ( bik0( r1 ) - bik1( r1 ) ) / hrad
        r1 = .2 * r1
        gmassh = gmassh +
     +         0.8 * r1**2 * ( bik0( r1 ) - bik1( r1 ) ) / ( 5. * hrad )
        gmassh = 2. * cmpmas( icmp ) * r * gmassh
c cubic halo - density tapers to zero at hrad
      else if( ctype( icmp ) .eq. 'CUBI' )then
        r1 = abs( r ) / hrad
        if( r1 .lt. 1.d0 )then
          gmassh = hmass * r1**3 * ( 5. * r1**3 - 9. * r1**2 + 5. )
        else
          gmassh = hmass
        end if
      else if( numcal .or. fixed )then
c numerical potentials or models with fixed rotation curves
        gmassh = vcirc( r )**2
        jcmp = icmp
        do icmp = 1, ncmp
          if( disc( icmp ) )gmassh = gmassh - vdisc( r )**2
        end do
        icmp = jcmp
        gmassh = max( r * gmassh, 0.d0 )
      end if
c none of the above - try numerical integral of halo density expression
      if( gmassh .lt. 0.d0 )then
        r1 = 0
        npts = 32
        allocate ( absc( npts ) )
        allocate ( wt( npts ) )
        call GLqtab( r1, r, npts, wt, absc )
        gmassh = 0
        do i = 1, npts
          r1 = absc( i )
          gmassh = gmassh + wt( i ) * gmfunh( r1 )
        end do
        gmassh = 4. * pi * gmassh
      end if
      return
      end
