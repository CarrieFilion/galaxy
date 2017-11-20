      real*8 function dfmtz( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Toomre-Zang DF for Mestel's V=const disk
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
      real*8 arg, q
c
      q = dfcns( 1, icmp )
      arg = -E * ( q + 1 ) + q * log( Lz )
      dfmtz = dfnorm( icmp ) * exp( arg )
      return
      end
