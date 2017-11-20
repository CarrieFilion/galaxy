      real function zthick( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the disk scale height in model units
c
c calling argument - also in model units
      real r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c tapered thickness
c      zthick = max( 1. - dble( r / rtrunc )**2, 0.d0 )
c      zthick = z0init( 1 ) * sqrt( zthick )
c constant scale height disc
      zthick = z0init( icmp )
      return
      end
