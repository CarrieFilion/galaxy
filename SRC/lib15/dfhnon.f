      real*8 function dfhnon( E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Henon's isotropic DF for the spherical isochrone (BT87 p239, eq 4-144)
c  apparently neither BT's expression nor Henon's integrates to the density!
c
c calling arguments
      real*8 E
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      real*8 re
      include 'inc/pi.f'
c
c the value of E tilde in BT is simply -E in my normalization
      re = sqrt( -E )
      dfhnon = 27. + 66. * E + 320. * E**2 + 240. * E**3 + 64. * E**4
     +         + 3. * ( 16. * E**2 - 28. * E - 9. ) *
     +                        asin( re ) / sqrt( E * ( E - 1. ) )
      dfhnon = re * dfhnon / ( 2. * ( E - 1. ) )**4
c normalization
      dfhnon = dfhnon / ( sqrt( 2.d0 ) * ( 2. * pi )**3 )
c$$$c Henon's expression (Ann d'Astrophys, v23, p474)
c$$$c is his E in my units?
c$$$      dfhnon = ( 4. * E**4 + 14. * E**3 + 14. * E**2 - 53. * E + 48. ) *
c$$$     +          sqrt( 1. - E ) / ( 1. + E )**4 +
c$$$     +   3. * ( 4. * E**2 - 22. * E + 9. ) * acos( E ) / ( 1. + E )**4.5
c$$$c normalization
c$$$      dfhnon = dfhnon / ( 4. * sqrt( 2.d0 ) * pi**3 )
      return
      end
