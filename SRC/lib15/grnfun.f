      real function grnfun( rs, rf, th, zf )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the softened gravitational potential at the given field point
c   due to a unit mass at the source point
c   The source point is assumed to lie on the line th = zs = 0
c
c calling arguments
      real rf, rs, th, zf
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      real sftpot
c
c local variable
      real d2
c
      if( p2d .or. p3d )then
        d2 = rs * rs + rf * rf - 2. * rs * rf * cos( th ) + zf * zf
        grnfun = sftpot( d2 )
      else
        call crash( 'GRNFUN', 'Unrecognized code' )
      end if
      return
      end
