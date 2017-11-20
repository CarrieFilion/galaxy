      subroutine wcell( x, y, z, w, rho )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Calculates the weights (in natural units) to be assigned to the cell
c   corners from a smooth density function.
c Note: the input coordinates are for the base corner of the cell in grid
c   units, while the external density function, rho, and the returned
c   weights are in natural units!
c
c calling arguments
      real x, y, z, w( * )
      real*8 rho
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c externals
      real uofgr, grofu
c
c local arrays
      integer n, mt
      parameter ( n = 4, mt = 10 )
      real*8 a( n ), c( mt ), q( 8 ), s( mt ), wg( mt ), wt( n )
c
c local variables
      integer i, j, k, l, m, nt
      logical first
      real*8 aa, afac, bb, bfac, d, dr, r, r1, r2, t, x2, y2, z2
      include 'inc/pi.f'
      save a, c, first, nt, s, wg, wt
c
      data first / .true. /
c
c compute factors for low order Gauss-Legendre quadrature once only
      if( first )then
        call GLqtab( 0.d0, 1.d0, n, wt, a )
c take account of cell size in natural units
        if( .not. s3d )then
          do i = 1, n
            wt( i ) = wt( i ) / lscale
          end do
        else
c weights and abscissae for polar angle integral for s3d
          nt = 1
          if( s3lmax .gt. 0 )then
            nt = mt
c assume mass is symmetric about z=0 and double weights
            d = .5 * pi
            call GLqtab( 0.d0, d, nt, wg, c )
            do i = 1, nt
              t = c( i )
              c( i ) = cos( t )
              s( i ) = sin( t )
              wg( i ) = 2. * wg( i )
            end do
          end if
        end if
        first = .false.
      end if
c
      if( c3d )then
c volume integrals of the density over a unit cube
        do i = 1, 8
          q( i ) = 0
        end do
        do k = 1, n
          z2 = ( z + a( k ) ) / lscale
          do j = 1, n
            y2 = ( y + a( j ) ) / lscale
            do i = 1, n
              x2 = ( x + a( i ) ) / lscale
              d = rho( x2, y2, z2 ) * wt( i ) * wt( j ) * wt( k )
              q( 1 ) = q( 1 ) + d
              q( 2 ) = q( 2 ) + a( i ) * d
              q( 3 ) = q( 3 ) + a( j ) * d
              q( 4 ) = q( 4 ) + a( k ) * d
              q( 5 ) = q( 5 ) + a( j ) * a( k ) * d
              q( 6 ) = q( 6 ) + a( i ) * a( k ) * d
              q( 7 ) = q( 7 ) + a( i ) * a( j ) * d
              q( 8 ) = q( 8 ) + a( i ) * a( j ) * a( k ) * d
            end do
          end do
        end do
c set weights - CIC scheme
c the following is the key to the various w's
c w( i ) = w(x,y,z) follows, with cell size of unit length
c
c w( 1 ) = w( 0,0,0 )           5 __________7
c w( 2 ) = w( 1,0,0 )            /.        /|
c w( 3 ) = w( 0,1,0 )           /_________/ |
c w( 4 ) = w( 1,1,0 )          6| :      8| |     z
c w( 5 ) = w( 0,0,1 )           | .       | |     |
c w( 5 ) = w( 0,0,1 )           | ........|.|3    |___y
c w( 6 ) = w( 1,0,1 )           | '1      | /     /
c w( 7 ) = w( 0,1,1 )           |'________|/     /x
c w( 8 ) = w( 1,1,1 )           2         4
c
        w( 1 ) = q( 1 ) - q( 2 ) - q( 3 ) - q( 4 ) +
     +           q( 5 ) + q( 6 ) + q( 7 ) - q( 8 )
        w( 2 ) = q( 2 ) - q( 7 ) - q( 6 ) + q( 8 )
        w( 3 ) = q( 3 ) - q( 5 ) - q( 7 ) + q( 8 )
        w( 4 ) = q( 7 ) - q( 8 )
        w( 5 ) = q( 4 ) - q( 6 ) - q( 5 ) + q( 8 )
        w( 6 ) = q( 6 ) - q( 8 )
        w( 7 ) = q( 5 ) - q( 8 )
        w( 8 ) = q( 8 )
c
      else if( p3d )then
c In cylindrical polar coordinates, the integration element is r d(phi) dr dz.
c   The phi integral is trivial because the mass distribution is axisymmetric
c   The arguments for the density function use the convention x=r, y=0.
        i = uofgr( x )
        dr = grofu( real( i + 1 ) ) - grofu( real( i ) )
c sum weighted contributions
        do i = 1, 4
          q( i ) = 0
        end do
        do k = 1, n
          z2 = ( z + a( k ) * dzg ) / lscale
          y2 = 0.d0
          do i = 1, n
            x2 = ( x + a( i ) * dr ) / lscale
            d = x2 * rho( x2, y2, z2 ) * wt( i ) * wt( k )
            q( 1 ) = q( 1 ) + d
            q( 2 ) = q( 2 ) + a( i ) * d
            q( 3 ) = q( 3 ) + a( k ) * d
            q( 4 ) = q( 4 ) + a( i ) * a( k ) * d
          end do
        end do
c renormalize for the cell dimensions - two lscale factors were included in
c   the wt factors above and alpha = 2 pi / na
        do i = 1, 4
          q( i ) = alpha * q( i ) * dzg * dr
        end do
c set weights - CIC scheme
        w( 1 ) = q( 1 ) - q( 2 ) - q( 3 ) + q( 4 )
        w( 2 ) = q( 2 ) - q( 4 )
        w( 5 ) = q( 3 ) - q( 4 )
        w( 6 ) = q( 4 )
c
      else if( p3a )then
c On the axisymetric polar grid, the volume integral is r d(phi) dr dz.
c   The phi integral is just 2pi because the grid is axisymmetric.
c   The arguments for the density function use the convention x=r, y=0.
        i = uofgr( x )
        dr = grofu( real( i + 1 ) ) - grofu( real( i ) )
c sum weighted contributions
        do i = 1, 4
          q( i ) = 0
        end do
        do k = 1, n
          z2 = ( z + a( k ) * dzg ) / lscale
          y2 = 0.d0
          do i = 1, n
            x2 = ( x + a( i ) * dr ) / lscale
            d = x2 * rho( x2, y2, z2 ) * wt( i ) * wt( k )
            q( 1 ) = q( 1 ) + d
            q( 2 ) = q( 2 ) + a( i ) * d
            q( 3 ) = q( 3 ) + a( k ) * d
            q( 4 ) = q( 4 ) + a( i ) * a( k ) * d
          end do
        end do
c renormalize for the cell dimensions - two lscale factors were included in
c   the wt factors above
        do i = 1, 4
          q( i ) = 2. * pi * q( i ) * dzg * dr
        end do
c set weights - CIC scheme
        w( 1 ) = q( 1 ) - q( 2 ) - q( 3 ) + q( 4 )
        w( 2 ) = q( 2 ) - q( 4 )
        w( 3 ) = q( 3 ) - q( 4 )
        w( 4 ) = q( 4 )
c
      else if( s3d )then
c Spherical grid: the volume element is sin(theta) d(theta) d(phi) r^2 dr
c   The phi integral just 2pi when the mass is axisymmetric and we need
c   evaluate only the m=0 terms
c   The calling arguments use the convention x=r, y=z=0
        do l = 0, s3lmax
          m = 2 * l
          w( m + 1 ) = 0
          w( m + 2 ) = 0
        end do
c calculate the actual radial extent of the cell with inner radius x
        i = uofgr( x ) + .5
        r1 = s3rad( i + 1 )
        r2 = s3rad( i + 2 )
        dr = r2 - r1
        y2 = 0
c integrate over radius
        do i = 1, n
          r = x + a( i ) * dr
c integrate over polar angle
          do j = 1, nt
            if( nt .gt. 1 )then
              x2 = r * s( j ) / lscale
              z2 = r * c( j ) / lscale
              d = 2. * pi * rho( x2, y2, z2 ) *
     +            wg( j ) * wt( i ) * dr * s( j )
c get Plms
              call s3tplm( c( j ), .false. )
            else
              x2 = r / lscale
              z2 = 0
              d = 4. * pi * rho( x2, y2, z2 ) *
     +            wt( i ) * dr
              call s3tplm( 0.d0, .false. )
            end if
c work over l - extra factor of r**2 for the radial interval
            aa = r**2 / r2
            afac = r / r2
            bb = r
            bfac = r1 / r
            k = 1
            do l = 0, s3lmax
              k = k + l
c only even l values give non-zero contributions
              if( mod( l, 2 ) .eq. 0 )then
                m = 2 * l
c interior mass
                w( m + 1 ) = w( m + 1 ) + d * aa * plm( k )
c exterior mass
                w( m + 2 ) = w( m + 2 ) + d * bb * plm( k )
              end if
              aa = aa * afac
              bb = bb * bfac
            end do
          end do
        end do
      else
        call crash( 'WCELL', 'Unrecognized grid' )
      end if
      return
      end
