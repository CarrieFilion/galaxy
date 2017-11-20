      real*8 function dfevco( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      real*8 E, Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      integer n, r, r2
      real*8 a, arg, b, dfm
      include 'inc/pi.f'
c
      dfm = dfcns( 1, icmp )
c equation (2.7) of Evans & Collett (MNRAS v264 p353) for F(E) only
      if( dfm .lt. 0. )then
        arg = .5 * E
c double value to make all stars direct - except for zero angular momentum
        dfevco = 1. / ( pi * ( exp( -arg ) - exp( arg ) ) )**2
        if( Lz .eq. 0. )dfevco = .5 * dfevco
      else
        n = nint( dfm )
        if( n .eq. 0 )then
c equation (2.14) of Evans & Collett
          dfevco = exp( -E ) * ( exp( -.5 * Lz**2 ) + 1. )
        else
c equation (2.11) of Evans & Collett
          dfevco = exp( -E - .5 * Lz**2 )
          if( Lz .ne. 0. )then
            a = 2**n
            do r = 0, n
              r2 = 2 * r
              if( r .gt. 0 )a = a * dble( n + 1 - r )
     +                                         / dble( r2 * ( r2 - 1 ) )
              arg = E * dble( n + r + 1 )
              dfevco = dfevco + a * real( ( n + r + 1 )**( r + 1 ) )
     +                          * Lz**r2 * exp( -arg )
            end do
          end if
        end if
c double value to make all stars direct
        dfevco = dfevco / ( 2. * pi**2 )
c taper out radial orbits at high energies if requested
        if( indexi( icmp ) .ne. 0 )then
          a = ( Emaxe - E ) / ( Emaxe - Emine )
          b = Lz / Lztmx( icmp )
          arg = 4. * ( a**2 + b**2 )
          if( arg .lt. .25 )then
            dfevco = 0
          else if( arg .lt. 1. )then
            arg = sqrt( arg )
            dfevco = dfevco * cos( pi * arg )**2
          end if
        end if
      end if
      return
      end
