      real*8 function eminx( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns smallest value of -e that will keep a star within rmax in
c   a hard Kuz'min/Toomre disc
c
c calling argument
      real*8 x
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c Kuz'min/Toomre disc
      if( ctype( icmp ) .eq. 'KT  ' )then
        eminx = phimax**2 - ( x / rmax )**2
        eminx = max( eminx, 0.d+0 )
        eminx = .5 * ( sqrt( eminx ) - phimax )
      else
        call crash( 'EMINX', 'Unknown disc' )
      end if
      return
      end
