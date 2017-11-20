      real*8 function dfevan( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Evans functions for power law discs (Evans & Read 1998, MNRAS, v300, p83)
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
      real*8 dfb, dfm
c
c dfm = \gamma in their notation and dfb = ( 1 + \gamma ) / \beta
      dfm = dfcns( 1, icmp )
      dfb = dfcns( 2, icmp )
      dfevan = dfnorm( icmp ) * ( Lz / sqrt( abs( E ) ) )**dfm
      dfevan = dfevan * abs( evbeta * E )**dfb
      return
      end
