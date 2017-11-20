      real*8 function sphrfn( v )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Density function which gives the weighted phase space density for SPHRHO.
c   Since the DF is assumed a function of E alone, a single integration only
c   is required and the weight factor is simply v^2.
c The external requires two arguments, but the angular momentum is ignored
c   for the isotropic DFs assumed here
c The routine works entirely in natural units.
c
c calling argument
      real*8 v
c
c common block
c
      include 'inc/moment.f'
c
c external
      real*8 distfn
c
c local variable
      real*8 E
c
      E = pot + .5 * v * v
      sphrfn = v * v * distfn( E, 0.d0 )
      return
      end
