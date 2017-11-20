      real function vrtfor( rs, rf, th, zf )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the cylindrical vertical component of the softened gravitational
c   force at the given field point due to a unit mass at the source point
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
      include 'inc/grids.f'
c
c external
      real softm
c
c local variables
      real d2, x
c
      if( p3d )then
        d2 = rs * rs + rf * rf - 2. * rs * rf * cos( th ) + zf * zf
        if( tsoft .eq. 1 )then
c Plummer softening
          vrtfor = -zf / ( d2 + softl2 )**1.5
        else
c cubic spline softening
          vrtfor = 0
          if( zf .ne. 0. )then
            x = sqrt( d2 ) / softl
            vrtfor = -zf * softm( x ) / d2**1.5
          end if
        end if
      else
        call crash( 'VRTFOR', 'Unrecognized code' )
      end if
      return
      end
