      real*8 function vdisc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns circular velocity due to disc component only
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
c      real sandvel
      real*8 bik0, bik1, confhg, frdisc, frdgau, frdisi
c
c local variables
      integer ipar
      real*8 a, arg, cc, scale, x, x1, xi2
      include 'inc/pi.f'
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'VDISC',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        scale = cmpmas( icmp ) / rscale( icmp )
        cc = dfcns( 3, icmp )
      else
        call crash( 'VDISC', 'Not a disc component' )
      end if
      x = r / rscale( icmp )
c no disc
      if( ( r .eq. 0. ) .or. ( .not. disc( icmp ) ) )then
        vdisc = 0.
c Kuz'min/Toomre disc
      else if( ctype( icmp ) .eq. 'KT  ' )then
        vdisc = sqrt( scale ) * x / ( 1. + x * x )**.75
c softened Kuz'min/Toomre disc
      else if( ctype( icmp ) .eq. 'SOFT' )then
        vdisc = sqrt( scale ) * x / ( rstar * rstar + x * x )**.75
c Maclaurin/Freeman/Kalnajs disc (aka Hunter model 1)
      else if( ( ctype( icmp ) .eq. 'MFK ' ) .or.
     +         ( ( ctype( icmp ) .eq. 'HUNT' ) .and.
     +           ( impar( icmp ) .eq. 1 ) ) )then
        if( abs( x ) .le. 1. )then
          vdisc = sqrt( scale * 3. * pi / 4. ) * x
        else
          a = 1. / abs( x )
          vdisc = 1.5 * x * x * ( asin( a ) - a * sqrt( 1. - a * a ) )
          vdisc = sqrt( scale * vdisc )
        end if
c radius is a dummy argument - Mestel/Toomre/Zang disc is self-similar
      else if( ctype( icmp ) .eq. 'MTZ ' )then
        vdisc = sqrt( scale )
c exponential disc
      else if( ctype( icmp ) .eq. 'EXP ' )then
        if( abs( x ) .lt. 1.e-10 )then
          vdisc = 3.416 * sqrt( scale ) * x
        else
          arg = .5 * abs( x )
          vdisc = .5 * ( bik0( arg ) - bik1( arg ) )
          if( abs( x ) .gt. 3. )vdisc = max( vdisc, 1. / abs( x )**3 )
          vdisc = x * sqrt( scale * vdisc )
        end if
c SC disc
      else if( ctype( icmp ) .eq. 'SC  ' )then
        vdisc = -frdisc( abs( x ) )
        vdisc = sqrt( vdisc * abs( x ) )
c isochrone disc
      else if( ctype( icmp ) .eq. 'ISOC' )then
        x1 = sqrt( 1. + x * x )
        vdisc = sqrt( scale ) * x / ( sqrt( x1 ) * ( 1. + x1 ) )
c Hunter models
      else if( ( ctype( icmp ) .eq. 'HUNT' ) .and.
     +         ( abs( x ) .lt. 1.d0 ) )then
        xi2 = 1. - x * x
        ipar = impar( icmp )
        if( ipar .eq. 2 )then
          vdisc = x * sqrt( 1.47262156 * scale * ( 1. + 3. * xi2 ) )
        else if( ipar .eq. 3 )then
          vdisc = x * sqrt( 1.28854386 * scale *
     +                                 ( 1. + 2. * xi2 + 5. * xi2**2 ) )
        else if( ipar .eq. 4 )then
          vdisc = x * sqrt( 1.20800987 * scale *
     +                   ( 1. + 1.8 * xi2 + 3. * xi2**2 + 7 * xi2**3 ) )
        else if( ipar .eq. 5 )then
          vdisc = x * sqrt( 1.1627095 * scale * ( 1. + 1.71428571 * xi2
     +             + 2.57142857 * xi2**2 + 4. * xi2**3 + 9. * xi2**4 ) )
        else
          call crash( 'VDISC', 'Hunter disc model not programmed' )
        end if
c Sanders modified gravity model with r0 = 8 * exp scale length
c      else if( ctype( icmp ) .eq. 10 )then
c        vdisc = sqrt( scale ) * sandvel( abs( x ) )
c power law disc
      else if( ctype( icmp ) .eq. 'POWE' )then
        vdisc = sqrt( scale ) * abs( x )**cc
c Gaussian disc - expression (32) from Toomre (1963, ApJ v138 p385)
      else if( ctype( icmp ) .eq. 'GAUS' )then
        if( epsiln .eq. 0.d0 )then
c   Note: Sigma(x) = M / ( 2 pi a^2 ) exp( - x^2 / 2a^2 )
          vdisc = scale * sqrt( .125 * pi ) *
     +                                 confhg( 1.5d0, 2.d0, -.5 * x**2 )
          vdisc = x * sqrt( vdisc )
        else
          vdisc = scale * frdgau( x )
          vdisc = sqrt( -x * vdisc )
        end if
c finite Mestel disc
      else if( ctype( icmp ) .eq. 'FMES' )then
        if( x .lt. 1. )then
          vdisc = sqrt( .5 * pi * scale )
        else
          vdisc = sqrt( scale * asin( 1. / x ) )
        end if
c Rybicki disc - equation (2.1) of Evans & Collett MNRAS v264 p353
      else if( ctype( icmp ) .eq. 'RYBI' )then
        vdisc = scale * ( 1. - 1. / sqrt( 1. + x * x ) )
        vdisc = sqrt( vdisc )
c Polynomial disc - equation (25) of Sawamura PASJ v40 p279
      else if( ctype( icmp ) .eq. 'SPLY' )then
        if( x .le. 1.d0 )then
          x1 = 1 - x * x
          vdisc = scale * 0.75 * pi * ( 1. + 2. * cc * x1 )
     +                                               / ( 1. + 0.4 * cc )
          vdisc = x * sqrt( vdisc )
        else
          vdisc = frdisc( r )
          vdisc = max( -x * vdisc, 0.d0 )
          vdisc = sqrt( vdisc )
        end if
c return default value for an unknown model
      else if( ctype( icmp ) .eq. 'UNKN' )then
        vdisc = sqrt( scale )
c composite model
      else if( ctype( icmp ) .eq. 'COMP' )then
        vdisc = -frdisc( r )
        vdisc = sqrt( x * vdisc )
c Donner-Thomasson double exponential disc  A&A v290 p785
      else if( ctype( icmp ) .eq. 'DOTH' )then
        if( abs( x ) .lt. 1.e-10 )then
          vdisc = 3.416 * sqrt( scale ) * x
        else
          arg = .5 * abs( x )
          vdisc = arg**2 * ( bik0( arg ) - bik1( arg ) )
          arg = 2 * arg
          vdisc = vdisc - .5 * arg**2 * ( bik0( arg ) - bik1( arg ) )
          vdisc = max( vdisc, 0.d0 )
          vdisc = sqrt( 8. * scale * vdisc / 3. )
        end if
      else
c none of the above - generic expression
        vdisc = frdisi( r )
        vdisc = max( -x * vdisc, 0.d0 )
        vdisc = sqrt( vdisc )
      end if
      vdisc = sign( vdisc, r )
      return
      end
