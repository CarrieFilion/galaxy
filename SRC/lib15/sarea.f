      real*8 function sarea( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the area enclosed in the meriodional plane by the zero velocity
c   curve for a given E & Lz in the model potential
c   It is needed as the Jacobian for the DF( E, Lz ) in a flattened spheroid
c   uses QUADPACK routine DQAGS
c
c calling agruments
      real*8 E, Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/orbval.f'
c
c externals
      external zmax
      real*8 quad_gnr, rcirc, zmax
c
c local variables
      integer ier
      real*8 eccen, epsa, epsr, r
      include 'inc/pi.f'
c
c get extreme radii in mid-plane
      call rlims( E, Lz )
c special case of particle at rest at centre
      eccen = 0
      if( rapo .gt. 1.d-10 )eccen = ( rapo - rperi ) / ( rapo + rperi )
      if( eccen .gt. .001 )then
c accurate area - symmetric about z = 0
        epsr = 1.d-5
        epsa = 1.d-5
        sarea = quad_gnr( zmax, rperi, rapo, epsa, epsr, ier )
        if( ier .ne. 0 )then
          print *, 'ier =', ier, ' from QUAD_GNR'
          call crash( 'SAREA', 'QUADPACK error' )
        end if
        sarea = 2. * sarea
      else if( eccen .gt. 0. )then
c assume an ellipse with semi-axes zmax( rcirc) and ( rapo - rperi ) / 2
        r = rcirc( Lz )
        sarea = pi * .5 * ( rapo - rperi ) * zmax( r )
      else
c circular orbit
        sarea = 0
      end if
      return
      end
