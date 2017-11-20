      real*8 function rhdsph( r )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c returns the mean spherical density at the given radius by integrating
c   over the polar angle
c   uses QUADPACK routine DQAGS
c
c calling argument in model units, not grid units
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      real*8 rs
      common / rhdspc / rs
c
c externals
      external rhdspi
      real*8 quad_gnr
c
c local variables
      integer ier
      real*8 epsa, epsr, zero
      include 'inc/pi.f'
      parameter ( zero = 0. )
c
c integrate over polar angle - half-the range only as we assume reflection symm
      rs = r
      epsa = 1.e-6
      epsr = 1.e-6
      rhdsph = quad_gnr( rhdspi, zero, .5 * pi, epsa, epsr, ier )
      if( ier .ne. 0 )then
        print *, 'ier =', ier, ' from QUAD_GNR'
        call crash( 'RHDSPH', 'QUADPACK error' )
      end if
      return
      end

      real*8 function rhdspi( theta )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c inner integrand over polar angle of the disk volume density at a point
c   specified in (r,theta) coordinates
c
      real*8 theta
c
c common block
c
      real*8 rs
      common / rhdspc / rs
c
c external
      real*8 rhorz
c
c local variables
      real*8 R, s, z
c
c convert spherical polars (r,theta) to cylindrical polars (R,z)
      s = sin( theta )
      z = rs * cos( theta )
      R = rs * s
c integrand needs a sin(theta) term
      rhdspi = s * rhorz( R, z )
      return
      end
