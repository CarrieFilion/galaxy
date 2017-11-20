      real function azavpt( r, z )
c  Copyright (C) 2015, Jerry Sellwood
c
c Returns azimuthally averaged potential on the grid.  Note the calling
c   arguments are in grid units but the value returned is in natural units!
c Similar function to the routine MEANPH, but the potential is defined by
c   POTGRD which computes the work done along some path
      use aarrays
      implicit none
c
c calling arguments - in mesh spaces
      real r, z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/grids.f'
c
c external
      real potgrd
c
c local variables
      integer i, j, k, m, n
      real f, rg, t, x, y
      include 'inc/pi.f'
      parameter ( m = 20, f = 2. * pi / m )
c
      if( r .lt. 0. )then
        print *, r, z
        call crash( 'AZAVPT', 'Impossible calling arguments' )
      end if
      if( analyt )call crash( 'AZAVPT', 'analyt true' )
      call switch( ngrid )
c the mesh values are the values of the potential from the Poisson solver
c   which is normalized to tend to zero at infinity.  potgrd is zero at the
c   origin.  ptoffs is a constant potential offset to make the two agree
c   at the grid edge
      if( iusepf .eq. 0 )then
        rg = rgrid( jgrid )
        if( c3d )then
c find average potential of four face centre values in the mid-plane
c   each one plane in from the boundary
          n = ngxy * ( ngz / 2 )
          i = n + ngx + ngx / 2 + 1
          j = n + ngx * ( ngy - 2 ) + ngx / 2 + 1
          k = n + ngx * ( ngy / 2 )
          ptoffs = .25 * ( grdfld( i, 4 ) - potgrd( 0., -ym, 0. )
     +                   + grdfld( j, 4 ) - potgrd( 0., ym, 0. )
     +               + grdfld( k + 2, 4 ) - potgrd( -xm, 0., 0. )
     +         + grdfld( k + ngx - 1, 4 ) - potgrd( xm, 0., 0. ) )
        else if( p3d )then
          n = nr( jgrid )
          n = ( ( n - 1 ) * ngz + ngz / 2 ) * na + 1
          ptoffs = grdfld( n, 4 ) - potgrd( rg, 0., 0. )
        else if( p3a )then
          n = ( ngz / 2 + 1 ) * nr( jgrid )
          ptoffs = grdfld( n, 4 ) - potgrd( rg, 0., 0. )
        else
c generic method
          ptoffs = -potgrd( rg, 0., 0. )
          oldc( 1, 1 ) = .9999 * rg
          oldc( 2, 1 ) = 0
          oldc( 3, 1 ) = 0
          iz( 1 ) = 1
          label( 1 ) = ngrid
          iflag( 1 ) = 2
          jgrid = ngrid
          call getacc( 1 )
          ptoffs = ptoffs + gpot( 1 )
        end if
        iusepf = 1
        ptoffs = ptoffs * gvfac**2
        print *, 'Potential offset =', ptoffs
      end if
c averaging needed only for c3d
      if( c3d )then
        if( r .eq. 0. )then
          azavpt = potgrd( 0., 0., z ) * gvfac**2 + ptoffs
        else
c average around a ring
          azavpt = 0
          do i = 1, m
            t = f * real( i - 1 )
            x = r * cos( t )
            y = r * sin( t )
            azavpt = azavpt + potgrd( x, y, z )
          end do
          azavpt = azavpt * gvfac**2 / real( m ) + ptoffs
        end if
      else
c mass distribution and forces are axisymmetric anyway
        azavpt = potgrd( r, 0., z ) * gvfac**2 + ptoffs
      end if
      return
      end
