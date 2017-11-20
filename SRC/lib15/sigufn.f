      real*8 function sigufn( xi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling argument
      real*8 xi
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gsigma, phitot, vcirc
c
      sigufn = xi**dfcns( 1, icmp ) * gsigma( xi ) *
     +        ( phitot( xi ) + vcirc( xi )**2 )
      return
      end
