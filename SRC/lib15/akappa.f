      real*8 function akappa( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns epicyclic frequency in mass model
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
      real*8 omegac, vcgrad
c
c local variables
      integer ip
      real*8 a, c, omega, r2
      include 'inc/pi.f'
c
c set nonsense result
      akappa = -1
c known expressions apply only for simple models
      if( ncmp .eq. 1 )then
        if( disc( 1 ) )then
c pure disc models
          if( ctype( 1 ) .eq. 'KT  ' )then
c Kuz'min/Toomre disc
            akappa = sqrt( 4. + r * r ) / ( 1. + r * r )**1.25
          else if( ctype( 1 ) .eq. 'SOFT' )then
c softened potential Kuz'min/Toomre disc
            akappa = sqrt( 4. * rstar * rstar + r * r ) /
     +                                   ( rstar * rstar + r * r )**1.25
          else if( ( ctype( 1 ) .eq. 'MFK ' ) .and.
     +             ( r .lt. 1.d0 ) )then
c Maclaurin/Freeman/Kalnajs disc - independent of radius within disc
            akappa = sqrt( 3. * pi )
          else if( ctype( 1 ) .eq. 'MTZ ' )then
c Mestel/Toomre/Zang disc
            a = .01
            a = max( a, r )
            akappa = sqrt( 2. ) / a
          else if( ctype( 1 ) .eq. 'ISOC' )then
c isochrone disc
            r2 = r * r
            a = sqrt( 1. + r2 )
            akappa = ( ( 4. + r2 ) * ( 1. + a ) + 2. * r2 ) /
     +                               ( ( 1. + r2 ) * a * ( 1. + a )**3 )
            akappa = sqrt( akappa )
          else if( ( ctype( 1 ) .eq. 'FMES' ) .and.
     +             ( r .lt. 1.d0 ) )then
c finite Mestel disc
            a = .01
            a = max( a, r )
            akappa = sqrt( 2. ) / a
          else if( ( ctype( 1 ) .eq. 'SPLY' ) .and.
     +             ( r .lt. 1.d0 ) )then
c Polynomial disc - equation (26) of Sawamura PASJ v40 p279
            c = dfcns( 3, 1 )
            akappa =  3. * pi * ( 1. + c * ( 2. - 3. * r * r ) )
     +                                                / ( 1. + 0.4 * c )
            akappa = sqrt( cmpmas( 1 ) * akappa )
          end if
        else
          if( ctype( 1 ) .eq. 'PLUM' )then
c Plummer model
            akappa = sqrt( 4. + r * r ) / ( 1. + r * r )**1.25
          else if( ctype( 1 ) .eq. 'JAFF' )then
c Jaffe model
            a = .00001
            a = max( a, r )
            akappa = sqrt( 2. + a ) / ( a * ( 1. + a ) )
          else if( ctype( 1 ) .eq. 'KEPL' )then
c Keplerian rotation curve
            akappa = r**( -1.5 )
          else if( ctype( 1 ) .eq. 'POWE' )then
c general power law rotation curve
            r2 = r
            c = dfcns( 3, 1 )
            if( c .lt. 1.d0 )r2 = max( r, 1.d-12 )
            akappa = sqrt( 2. * ( 1. + c ) ) * r2**( c - 1. )
          else if( ctype( 1 ) .eq. 'HERN' )then
c Hernquist halo
            akappa = ( 3. + r ) / ( 1. + r )**3
            a = 1.e-12
            a = max( a, r )
            akappa = sqrt( akappa / a )
          else if( ctype( 1 ) .eq. 'NFW ' )then
c NFW halo
            a = 1.e-5
            a = max( a, r )
            c = a / ( 1. + a )
            akappa = log( 1. + a ) - c + c**2
            akappa = sqrt( akappa / a**3 )
c uniform sphere
          else if( ctype( 1 ) .eq. 'UNIS' ) then
            if( r .gt. rscale( 1 ) )then
              akappa = cmpmas( 1 ) / r**3
            else
              akappa = 4. * cmpmas( 1 ) / rscale( 1 )**3
            end if
            akappa = sqrt( akappa )
          end if
        end if
      else if( fixed )then
        do ip = 1, ncmp
          if( ctype( ip ) .eq. 'FE  ' )then
c Fall/Efstathiou halo
            r2 = r / rscale( ip )
            akappa = dfcns( 3, ip ) * sqrt( 4. + 2. * r2 * r2 ) /
     +                               ( rscale( ip ) * ( 1. + r2 * r2 ) )
          end if
        end do
      end if
c none of the above - evaluate defining expression
      if( akappa .lt. 0.d0 )then
        omega = omegac( r )
        akappa = 2. * omega * ( omega + vcgrad( r ) )
        akappa = max( akappa, 0.d0 )
        akappa = sqrt( akappa )
      end if
      return
      end
