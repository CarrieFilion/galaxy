      real*8 function gsigma( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns value of surface density of disc component at the requested radius
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
c local variables
      integer ipar
      real*8 a, c, scale, x, x1
      include 'inc/pi.f'
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'GSIGMA',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        scale = cmpmas( icmp ) / rscale( icmp )**2
        x = r / rscale( icmp )
        c = dfcns( 3, icmp )
      else
        call crash( 'GSIGMA', 'Called for a non-disc component' )
      end if
c ensure radius is within the active region
      gsigma = 0
      if( ( r .gt. rhole ) .or. ( r .lt. rmax ) )then
c Kuz'min/Toomre disc
        if( ( ctype( icmp ) .eq. 'KT  ' ) .or.
     +      ( ctype( icmp ) .eq. 'SOFT' ) )then
          gsigma = scale / ( 2. * pi * ( 1. + x * x )**1.5 )
c Maclaurin/Freeman/Kalnajs disc
        else if( ctype( icmp ) .eq. 'MFK ' )then
          if( x .lt. 1.d0 )gsigma =
     +                 scale * ( 3. / ( 2. * pi ) ) * sqrt( 1. - x * x )
c Mestel/Toomre/Zang disc
        else if( ctype( icmp ) .eq. 'MTZ ' )then
          x1 = 1.d-5
          x1 = max( x1, x )
          gsigma = scale / ( 2. * pi * x1 )
c exponential disc
        else if( ( ctype( icmp ) .eq. 'EXP ' ) .or.
     +           ( ctype( icmp ) .eq. 'SAND' ) )then
          gsigma = scale * exp( -x ) / ( 2. * pi )
c SC disc - 1.5625 factor comes from .1 / ( lscale * ts**2 )
        else if( ctype( icmp ) .eq. 'SC  ' )then
          gsigma = 1.5625 * scale /
     +                    ( ( 1. + 20. * x )**.75 * ( 1. + .2 * x )**2 )
c isochrone disc
        else if( ctype( icmp ) .eq. 'ISOC' )then
          if( x .eq. 0.d0 )then
            gsigma = scale / ( 6. * pi )
          else
            a = sqrt( 1. + x * x )
            gsigma =
     +             scale * ( log( x + a ) - x / a ) / ( 2. * pi * x**3 )
          end if
c Hunter discs
        else if( ctype( icmp ) .eq. 'HUNT' )then
          ipar = impar( icmp )
          if( x .lt. 1.d0 )gsigma = scale * dble( 2 * ipar + 1 )
     +            * ( 1. - x * x )**( dble( ipar ) - .5 ) / ( 2. * pi )
c power law disc
        else if( ctype( icmp ) .eq. 'POWE' )then
          x1 = max( x, 1.d-10 )
          gsigma = scale * evsigma * x1**( 2. * c - 1. )
c "cosine" disc - c = r / a
        else if( ctype( icmp ) .eq. 'COSI' )then
          x1 = ( 1. - x ) * c
          if( abs( x1 ) .lt. 1. )gsigma =
     +                          scale * c * sqrt( 1. - x1 * x1 ) / pi**2
c Gaussian disc
        else if( ctype( icmp ) .eq. 'GAUS' )then
          gsigma = scale * exp( -.5 * x * x ) / ( 2. * pi )
c finite Mestel disc
        else if( ctype( icmp ) .eq. 'FMES' )then
          if( x .lt. 1.d0 )gsigma =
     +                 scale * ( 1. - 2. * asin( x ) / pi ) / ( 4. * x )
c arbitrary model
        else if( ctype( icmp ) .eq. 'ARBT' )then
          print *, 'Library routine called for arbitrary model'
c Rybicki disc - equation (2.2) of Evans & Collett MNRAS v264 p353
        else if( ctype( icmp ) .eq. 'RYBI' )then
          gsigma = scale / ( 2. * pi * sqrt( 1. + x * x ) )
c Polynomial disc - equation (24) of Sawamura PASJ v40 p279
        else if( ctype( icmp ) .eq. 'SPLY' )then
          if( x .lt. 1.d0 )then
            x1 = 1. - x * x
            gsigma = scale * sqrt( x1 ) *
     +   ( 9. + c * ( 16. * x1 - 6. ) ) / ( 6. * pi * ( 1. + 0.4 * c ) )
          end if
c composite disc
        else if( ctype( icmp ) .eq. 'COMP' )then
          x1 = 3 * x
          gsigma = scale *
     +      ( exp( -x ) / ( 2. * pi ) + .1 * 9. * exp( -x1 * x1 ) / pi )
c Donner-Thomasson double exponential model A&A v290 p785
        else if( ctype( icmp ) .eq. 'DOTH' )then
          gsigma = scale *
     +                 ( exp( -x ) - exp( -2. * x ) ) * 2. / ( 3. * pi )
        else
          call crash( 'GSIGMA', 'Unknown type of disc' )
        end if
      end if
      return
      end
