      real function azifor( rs, rf, th, zf )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the azimuthal component of the softened gravitational force
c   at the given field point due to a unit mass at the source point
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
      if( p2d .or. p3d )then
        d2 = rs * rs + rf * rf - 2. * rs * rf * cos( th ) + zf * zf
        if( tsoft .eq. 1 )then
c Plummer softening
          azifor = -rs * sin( th ) / ( d2 + softl2 )**1.5
        else
c cubic spline softening
          azifor = 0
          if( d2 .gt. 0. )then
            x = sqrt( d2 ) / softl
            azifor = -rs * sin( th ) * softm( x ) / d2**1.5
          end if
        end if
      else
        call crash( 'AZIFOR', 'Unrecognized grid' )
      end if
      return
      end
