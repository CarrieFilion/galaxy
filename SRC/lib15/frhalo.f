      real*8 function frhalo( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns radial force from halo mass
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
      real*8 frdisc, frtot, gmassh, bik0, bik1
c
c local variables
      integer jcmp
      real*8 a, cc, hmass, hrad, x
c
      hmass = cmpmas( icmp )
      hrad = rscale( icmp )
      cc = dfcns( 3, icmp )
c
      frhalo = 0
      if( ( r .gt. 0. ) .and. ( .not. disc( icmp ) ) )then
        x = r / hrad
        if( fixed )then
c models with fixed rotation curves
          jcmp = icmp
          frhalo = frtot( r )
          do icmp = 1, ncmp
            if( disc( icmp ) )frhalo = frhalo - frdisc( r )
          end do
          icmp = jcmp
c frozen disc model
        else if( ctype( icmp ) .eq. 'ASDI' )then
          jcmp = icmp
          do icmp = 1, ncmp
            if( disc( icmp ) .and. ( cmpmas( icmp ) .gt. 0. ) )frhalo =
     +                     cmpmas( jcmp ) * frdisc( r ) / cmpmas( icmp )
          end do
          icmp = jcmp
c uniform sphere
        else if( ctype( icmp ) .eq. 'UNIS' )then
          if( x .le. 1. )then
            frhalo = -hmass * x / hrad**2
          else
            frhalo = -hmass / r**2
          end if
c Plummer sphere
        else if( ctype( icmp ) .eq. 'PLUM' )then
          frhalo = -hmass * r / ( hrad * hrad + r * r )**1.5
c cored isothermal model
        else if( ctype( icmp ) .eq. 'ISOT' )then
          frhalo = -hmass * x / ( r * r + hrad * hrad )
c Jaffe model
        else if( ctype( icmp ) .eq. 'JAFF' )then
          frhalo = -hmass / ( r * ( hrad + r ) )
c Kepler
        else if( ctype( icmp ) .eq. 'KEPL' )then
          frhalo = -hmass / ( r * r )
c TEDS model
        else if( ctype( icmp ) .eq. 'TEDS' )then
          a = cc * hrad
          frhalo = .5 * r / ( hrad * hrad + r * r )**1.5 +
     +             .09 * r / ( a * a + r * r )**1.5
          frhalo = -hmass * frhalo / .59
c Kuz'min-Kutuzov oblate spheroid - value in mid-plane
        else if( ctype( icmp ) .eq. 'KKSP' )then
          x = sqrt( hrad * hrad + r * r )
          frhalo = -hmass * r / ( x * ( cc * hrad + x )**2 )
c Hernquist model
        else if( ctype( icmp ) .eq. 'HERN' )then
          frhalo = -hmass / ( hrad + r )**2
c use Plummer sphere expression as default for an unknown model
        else if( ctype( icmp ) .eq. 'UNKN' )then
          frhalo = -hmass * r / ( hrad * hrad + r * r )**1.5
c generalized Aguilar-Merritt model
        else if( ctype( icmp ) .eq. 'AGME' )then
          if( x .lt. 1. )then
            frhalo = -hmass * r**( 1. + cc ) / hrad**( 3. + cc )
          else
            frhalo = -hmass / r**2
          end if
c spherical isochrone
        else if( ctype( icmp ) .eq. 'ISOC' )then
          a = sqrt( hrad * hrad + r * r )
          frhalo = -hmass * r / ( a * ( hrad + a )**2 )
c NFW profile
        else if( ctype( icmp ) .eq. 'NFW ' )then
          a = r / hrad + 1.d-10
          frhalo = ( 1. / ( a * ( 1. + a ) ) - log( 1. + a ) / a**2 )
     +                                                 * hmass / hrad**2
c singular isothermal sphere
        else if( ctype( icmp ) .eq. 'SISP' )then
          a = r + 1.d-10
          frhalo = -hmass / a
c Miyamoto-Nagai spheroid - value in mid-plane
        else if( ctype( icmp ) .eq. 'MYNG' )then
          x = R**2 + ( 1. + cc )**2
          frhalo = -hmass * r / x**1.5
c halo of Donner-Thomasson A&A v290 p785 - potential of a double exp disk!
        else if( ctype( icmp ) .eq. 'EXPH' )then
          a = min( r, 10.d0 )
          a = .5 * ( a / hrad )
          frhalo = 0.2 * a**2 * ( bik0( a ) - bik1( a ) ) / hrad
          a = .2 * a
          frhalo = frhalo +
     +            0.8 * a**2 * ( bik0( a ) - bik1( a ) ) / ( 5. * hrad )
          if( r .gt. 10. )frhalo = frhalo * ( 10. / r )
          frhalo = -2. * cmpmas( icmp ) * frhalo / r
c cubic halo - density tapers to zero at hrad
        else if( ctype( icmp ) .eq. 'CUBI' )then
          x = abs( r ) / hrad
          if( x .lt. 1.d0 )then
            frhalo = -hmass *
     +                      ( 5. * x**4 - 9. * x**3 + 5. * x ) / hrad**2
          else
            frhalo = -hmass / r**2
          end if
        else
c generic model
          frhalo = -gmassh( r ) / r**2
        end if
      end if
      return
      end
