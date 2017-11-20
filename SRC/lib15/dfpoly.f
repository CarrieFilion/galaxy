      real*8 function dfpoly( e )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c DF for a spherical polytrope
c
c calling argument
      real*8 e
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variable
      real*8 Phi1
c
      Phi1 = -cmpmas( icmp ) / rscale( icmp )
      if( e .lt. Phi1 )then
        dfpoly = dfnorm( icmp ) * ( -e + Phi1 )**dfcns( 1, icmp )
      else
        dfpoly = 0
      end if
      return
      end
