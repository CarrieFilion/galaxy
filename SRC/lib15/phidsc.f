      real*8 function phidsc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c value of potential in the disc plane from disc component only
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
c      real sandpot
      real*8 phgend, phiexp, phigau, phimes, phisc
c phdsci,
c
c local variables
      real*8 cc, scale, x, x1
      include 'inc/pi.f'
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'PHIDSC',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        scale = cmpmas( icmp ) / rscale( icmp )
        cc = dfcns( 3, icmp )
      else
        call crash( 'PHIDSC', 'Called for a non-disc component' )
      end if
      x = r / rscale( icmp )
c Kuz'min/Toomre disc
      if( ctype( icmp ) .eq. 'KT  ' )then
        phidsc = -scale / sqrt( 1. + x * x )
c softened Kuz'min/Toomre disc
      else if( ctype( icmp ) .eq. 'SOFT' )then
        phidsc = -scale / sqrt( rstar * rstar + x * x )
c Maclaurin/Freeman/Kalnajs disc
      else if( ctype( icmp ) .eq. 'MFK ' )then
        if( x .le. 1. )then
          phidsc = -scale * ( 3. * pi / 4. ) * ( 1. - .5 * x * x )
        else
          phidsc = -.75 * scale * ( ( 2. - x * x ) * asin( 1. / x )
     +                              + sqrt( x * x - 1. ) )
        end if
c Mestel/Toomre/Zang disc
      else if( ctype( icmp ) .eq. 'MTZ ' )then
        x1 = 1.e-20
        x1 = max( x, x1 )
        phidsc = scale * log( x1 )
c infinite exponential disc
      else if( ctype( icmp ) .eq. 'EXP ' )then
        phidsc = scale * phiexp( x )
c SC disc
      else if( ctype( icmp ) .eq. 'SC  ' )then
        phidsc = scale * phisc( x )
c isochrone disc
      else if( ctype( icmp ) .eq. 'ISOC' )then
        phidsc = -scale / ( 1. + sqrt( 1. + x * x ) )
c Sanders modified gravity model with r0 = 8*exp scale length
c      else if( ctype( icmp ) .eq. 'SAND' )then
c        phidsc = scale * sandpot( x )
c power law disc
      else if( ctype( icmp ) .eq. 'POWE' )then
        if( cc .lt. 0.d0 )then
          x1 = 1.e-20
          x1 = max( x, x1 )
        else
          x1 = x
        end if
        phidsc = scale * x1**( 2. * cc ) / ( 2. * cc )
c Gaussian disc
      else if( ctype( icmp ) .eq. 'GAUS' )then
        phidsc = scale * phigau( x )
c finite Mestel disc
      else if( ctype( icmp ) .eq. 'FMES' )then
        phidsc = scale * phimes( x )
c Rybicki disc - equation (2.2) of Evans & Collett MNRAS v264 p353
c   which is in fact the full 3-D expression
      else if( ctype( icmp ) .eq. 'RYBI' )then
        phidsc = scale * log( 1. + sqrt( 1. + x * x ) )
c Polynomial disc - equation (22) of Sawamura PASJ v40 p279
      else if( ctype( icmp ) .eq. 'SPLY' )then
        if( r .le. 1.d0 )then
          x1 = 1 - x * x
          phidsc = -scale * 0.375 * pi
     +                     * ( 1. + x1 + cc * ( x1 * x1 + 1. / 3. ) )
     +                     / ( 1. + 0.4 * cc )
        else
c          print *, 'r > 1 in polynomial disc - error in PHIDSC'
c          stop
          phidsc = -scale / x
        end if
c Donner-Thomasson double exponential disc A&A v290 p785
      else if( ctype( icmp ) .eq. 'DOTH' )then
c this expression does not give reasonable results!  Don't know what is wrong.
        phidsc = 4. * scale *
     +                      ( phiexp( x ) - 2. * phiexp( 2. * x ) ) / 3.
      else
c none of the above - evaluate potential directly from the surface density
        phidsc = phgend( x )
      end if
      return
      end
