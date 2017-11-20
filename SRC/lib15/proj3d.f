      real function proj3d( pos, idir )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns projected density of a mode from an expansion in spherical Bessel
c   functions - part of mode fitting software
c
c The projection direction is controlled by the calling argument idir and
c   the remaining two coordinates are held fixed at the values set in the
c   calling argument pos
c
c calling arguments
      integer idir
      real pos( 3 )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/anlys.f'
c
      common / sphlmn / l, m, n, mode
      integer l, m, n, mode
c
c external
      complex sphfn
c
c local arrays
      integer npts
      parameter ( npts = 10 )
      real*8 absc( npts ), wt( npts )
c
c local variables
      complex f
      integer i, ip
      real x, y, z
      real*8 a, b
c
      if( mod( jp, jnmax ) .ne. 1 )call crash( 'PROJ3D',
     +                     'First function is not lowest radial order' )
c compute integration limits
      if( idir .eq. 1 )then
        x = pos( 3 )
        y = pos( 2 )
        b = sqrt( 1. - x**2 - y**2 )
      end if
      if( idir .eq. 2 )then
        x = pos( 3 )
        z = pos( 1 )
        b = sqrt( 1. - x**2 - z**2 )
      end if
      if( idir .eq. 3 )then
        y = pos( 2 )
        z = pos( 1 )
        b = sqrt( 1. - y**2 - z**2 )
      end if
      proj3d = 0
      if( b .gt. 0.d0 )then
        a = -b
c compute weights and abscissae - Gauss-Legendre quadrature
        call GLqtab( a, b, npts, wt, absc )
c evaluate integral
        do i = 1, npts
          if( idir .eq. 1 )z = absc( i )
          if( idir .eq. 2 )y = absc( i )
          if( idir .eq. 3 )x = absc( i )
c sum over functions
          n = 0
          l = jp / jnmax + m
          do ip = jp, kp
            n = n + 1
            if( n .gt. jnmax )then
              n = 1
              l = l + 2
            end if
            f = sphfn( x, y, z )
            proj3d = proj3d + ( alp( ip, mode ) * real( f ) +
     +                          bet( ip, mode ) * aimag( f ) ) * wt( i )
          end do
        end do
        proj3d = proj3d / rbess**3
      end if
      return
      end
