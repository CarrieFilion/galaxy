      real*8 function dfomeg( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Kalnajs's DF for the Maclaurin disk
c
c calling arguments
      real*8 E, lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      dfomeg = dfnorm( icmp ) /
     +                  sqrt( omeg21 - 2. * ( E - Phi0 - omegak * Lz ) )
      return
      end
