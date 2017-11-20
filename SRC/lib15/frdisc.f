      real*8 function frdisc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c radial force from active disc
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
      real*8 bik0, bik1, confhg, frdgau, frdisi
c
c local variables
      real*8 a, arg, cc, scale, x, x1
      include 'inc/pi.f'
c
      scale = cmpmas( icmp ) / rscale( icmp )**2
      cc = dfcns( 3, icmp )
      x = r / rscale( icmp )
c potential well must be symmetric about r = 0
      if( ( r .eq. 0. ) .or. ( .not. disc( icmp ) ) )then
        frdisc = 0.
c Kuz'min/Toomre disc
      else if( ctype( icmp ) .eq. 'KT  ' )then
        frdisc = -scale * x / ( 1. + x * x )**1.5
c softened potential Kuz'min/Toomre disc
      else if( ctype( icmp ) .eq. 'SOFT' )then
        frdisc = -scale * x / ( rstar * rstar + x * x )**1.5
c Maclaurin/Freeman/Kalnajs disc
      else if( ctype( icmp ) .eq. 'MFK ' )then
        if( x .le. 1.d0 )then
          frdisc = -scale * ( 3. * pi / 4. ) * x
        else
          a = 1. / x
          frdisc = -scale * 1.5 * x *
     +                            ( asin( a ) - a * sqrt( 1. - a * a ) )
        end if
c Mestel/Toomre/Zang disc
      else if( ctype( icmp ) .eq. 'MTZ ' )then
        frdisc = -scale / x
c exponential disc
      else if( ctype( icmp ) .eq. 'EXP ' )then
        arg = .5 * abs( x )
        frdisc = -scale * arg * ( bik0( arg ) - bik1( arg ) )
c force of SC disc in 'natural' units
      else if( ctype( icmp ) .eq. 'SC  ' )then
        frdisc = -scale * 3.825 * x
     +               / ( ( 1. + 2. * x )**1.75 * ( 1. + .2 * x )**1.25 )
c isochrone disc
      else if( ctype( icmp ) .eq. 'ISOC' )then
        x1 = sqrt( 1. + x * x )
        frdisc = -scale * x / ( x1 * ( 1. + x1 )**2 )
c      end if
c Sanders modified gravity model with r0 = 8 * exp scale length
c      if( ctype( icmp ) .eq. 'SAND' )then
c        frdisc = -scale * sandvel( x )**2 / x
c power law disc
      else if( ctype( icmp ) .eq. 'POWE' )then
        frdisc = -scale * x**( 2. * cc - 1. )
c Gaussian disc
      else if( ctype( icmp ) .eq. 'GAUS' )then
        if( epsiln .eq. 0.d0 )then
          frdisc = -scale * sqrt( .125 * pi ) * x *
     +                                 confhg( 1.5d0, 2.d0, -.5 * x**2 )
        else
          frdisc = scale * frdgau( x )
        end if
c finite Mestel disc
      else if( ctype( icmp ) .eq. 'FMES' )then
        if( x .lt. 1.d0 )then
          frdisc = -scale * .5 * pi / x
        else
          frdisc = -scale * asin( 1. / x ) / x
        end if
c Rybicki disc - equation (2.1) of Evans & Collett MNRAS v264 p353
      else if( ctype( icmp ) .eq. 'RYBI' )then
        frdisc = -scale * ( 1. - 1. / sqrt( 1. + x * x ) ) / x
c Polynomial disc - equation (25) of Sawamura PASJ v40 p279
      else if( ( ctype( icmp ) .eq. 'SPLY' ) .and. ( x .le. 1.d0 ) )then
        x1 = 1 - x * x
        frdisc = 0.75 * pi * ( 1. + 2. * cc * x1 ) / ( 1. + 0.4 * cc )
        frdisc = -scale * x * frdisc
c composite model
      else if( ctype( icmp ) .eq. 'COMP' )then
c exponential part
        if( abs( x ) .lt. 1.e-10 )then
          frdisc = 3.416**2
        else
          arg = .5 * abs( x )
          frdisc = .5 * ( bik0( arg ) - bik1( arg ) )
          if( abs( x ) .gt. 3. )frdisc = max( frdisc, 1. / abs( x )**3 )
        end if
c Gaussian part - expression (32) from Toomre (1963, ApJ v138 p385)
c   Note: Sigma(0) = 1 / pi
        x1 = 3 * x
        frdisc = frdisc +
     +                  2.7 * sqrt( pi ) * confhg( 1.5d0, 2.d0, -x1**2 )
        frdisc = -x * scale * frdisc
c Donner-Thomasson exponential disc A&A v290 p785
      else if( ctype( icmp ) .eq. 'DOTH' )then
        arg = .5 * abs( x )
        frdisc = arg**2 * ( bik0( arg ) - bik1( arg ) )
        arg = 2 * arg
        frdisc = frdisc - .5 * arg**2 * ( bik0( arg ) - bik1( arg ) )
        frdisc = -8. * scale * frdisc / ( 3. * r )
c none of the above - generic expression
      else
        frdisc = frdisi( r )
c        call crash( 'FRDISC', 'Unknown type of disc' )
      end if
      return
      end
