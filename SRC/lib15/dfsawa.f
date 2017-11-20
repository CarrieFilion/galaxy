      real*8 function dfsawa( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Sawamura's DFs for his polynomial discs - PASJ v40 p 279
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
      integer n
      real*8 a, c, g, Phi, xi
      include 'inc/pi.f'
c
c r determined by \Phi(r) = E - .5 Lz^2 - inversion of equation (22)
      Phi = ( E - .5 * Lz**2 )
      g = dfcns( 3, icmp )
      if( g .eq. 0.d0 )then
        xi = 1. + Phi / ( .375 * pi )
      else
        c = 1. + g / 3. + Phi * ( 1. + .4 * g ) / ( .375 * pi )
        xi = ( -1. + sqrt( 1. - 4. * g * c ) ) / ( 2. * g )
      end if
c
      n = dfcns( 1, icmp ) + .1
c equation (27)
      if( n .eq. 0 )then
        dfsawa = 2. * ( 9. - 6. * g + 32. * g * xi )
     +                           / ( 9. * pi**3 * ( 1. + 2. * g * xi ) )
      else if( n .eq. 3 )then
c equation (28)
        c = g * xi
        a = 8. * ( 1. + .4 * g ) * Lz**2 / ( pi * ( 1. + 2. * c )**2 )
        dfsawa = xi**3 * ( 9. - 6. * g + 20. * c )
     +         + xi**2 * a * ( 27. - 18. * g + ( 116. - 24. * g ) * c
     +                          + 120. * c**2 )
     +         + 2. * xi * a**2 * ( 9. - 6. * g + ( 49. - 6. * g ) * c
     +                          + 80. * c**2 + 40. * c**3 ) / 3.
     +         + 2. * a**3 * ( 9. - 6. * g + ( 26. + 36. * g ) * c
     +                          - ( 54. - 36. * g ) * c**2 - 160. * c**3
     +                          - 80. * c**4 ) / 135.
        dfsawa = dfsawa * 8. / ( 9. * pi**3 * ( 1. + 2. * c ) )
      else
        call crash( 'DFSAWA', 'Unknown expression for DF const' )
      end if
      return
      end
