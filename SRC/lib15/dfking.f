      real*8 function dfking( e )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c DF for an isotropic King model
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
      real*8 eps
c
      eps = -e - cmpmas( icmp ) / rtidal
      if( eps .gt. 0.d0 )then
        dfking = dfnorm( icmp ) * ( exp( eps / sigmak**2 ) - 1.d0 )
      else
        dfking = 0
      end if
      return
      end
