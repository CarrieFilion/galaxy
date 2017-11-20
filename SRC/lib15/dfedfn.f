      real*8 function dfedfn( Psi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for Eddington's inversion formula - called from DFEDDI & DFNFW
c    see notes in  /home/sellwood/docs/models/Eddinv.tex
c
c calling argument
      real*8 Psi
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 frhalo, gmassh, rhohal, rofPsi
c
c local variables
      real*8 d1, d2, gm, haloc, hrad, r, rho, term, tiny
      real*8 x, x1, xm
      include 'inc/pi.f'
      parameter ( tiny = 1.d-26 )
c
c part that is independent of the model
      r = rofPsi( Psi )
      hrad = rscale( icmp )
      x = r / hrad
      if( ctype( icmp ) .eq. 'NFW ' )then
        if( x .gt. 1.1d-3 )then
c see notes
          x1 = 1. + x
          term = 2. * x1**2 * log( x1 ) - x * ( 2. + 3. * x )
          term = term / ( x1**2 * log( x1 ) - x * x1 )
          d1 = -( 1. + 4. * x + 3. * x**2 ) / x1**2
          d2 = 2. *
     +     ( 1. + 6. * x + 15. * x**2 + 16. * x**3 + 6. * x**4 ) / x1**4
          dfedfn = x * ( d2 + term * d1 ) / ( x1 * log( x1 ) - x )**2
        else
c approximate behavior for very small arguments
          dfedfn = 8.018 / ( x**3 + tiny )
        end if
        dfedfn = dfedfn / ( 4. * pi * hrad * cmpmas( icmp ) )
      else if( ctype( icmp ) .eq. 'EINA' )then
        gm = gmassh( r )
        rho = rhohal( r )
        term = 2. - 4. * pi * r**3 * rho / ( gm + tiny )
c Einasto halo
        xm = x + tiny
        haloc = dfcns( 3, icmp )
        d1 = -2.d0 * rho * x**haloc / ( hrad * xm )
        d2 = -d1 * ( 2.d0 * x**haloc - haloc + 1.d0 ) / ( hrad * xm )
        term = max( d2 + d1 * term / ( r + tiny ), 0.d0 )
        dfedfn = term / ( frhalo( r )**2 + tiny )
      else if( ctype( icmp ) .eq. 'ISOT' )then
c Isothermal halo - integrand can be written in closed form
        dfedfn = ( 9. + x**2 ) /
     +                 ( pi * cmpmas( icmp ) * hrad * ( 1. + x**2 )**2 )
      else if( ctype( icmp ) .eq. 'BURK' )then
c Burkert halo - Eddington inversion fails as the DF is -ve for this halo
        if( x .gt. 1.d-3 )then
          xm = ( 1. + x ) * ( 1. + x**2 )
          rho = rhohal( 0.d0 )
          d1 = -rho * ( 1. + 2. * x + 3. * x**2 ) / ( hrad * xm**2 )
          d2 = -rho * ( xm**2 * ( 2. + 6. * x ) -
     +            ( 1. + 2. * x + 3. * x**2 )**2 ) / ( hrad**2 * xm**4 )
          gm = gmassh( r )
          rho = rhohal( r )
          term = 2. / ( r + tiny ) - 4.* pi * r**2 * rho / ( gm + tiny )
          dfedfn = ( d2 + d1 * term ) / ( frhalo( r )**2 + tiny )
        else
c approximate behavior for very small arguments
          dfedfn = 0.7162 / ( x**3 + tiny )
        end if
      else
        call crash( 'DFEDFN', 'Unrecognized halo' )
      end if
      return
      end
