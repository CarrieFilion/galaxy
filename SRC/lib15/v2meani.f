      real*8 function v2meani( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes v2mean from a double integral over distribution fn
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
c externals
      real*8 velint
c
c set power factors for velocity weighting
      iu = 0
      iv = 2
      v2meani = velint( r )
      if( v2meani .ne. 0.d0 )then
c normalise by active surface density
        iu = 0
        iv = 0
        v2meani = v2meani / velint( r )
      end if
c in spherical models this integral includes both the theta and phi components
      if( ( .not. disc( icmp ) ) .and.
     +    ( .not. sphrod( icmp ) ) )v2meani = .5d0 * v2meani
      return
      end
