      subroutine massdf( rho )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to assign a mass distribution to the grid from a smooth density
c   function of position.  The external density function can be either
c   analytic or defined by a numerical integral of a DF over the velocities
      use aarrays
      implicit none
c
c calling argument
      external rho
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c external
      real grofu
c
c local array
      integer lmax, nw
c      parameter ( lmax = 16, nw = max( 2 * ( lmax + 1 ), 8 ) )
      parameter ( lmax = 16, nw = 34 )
      real wt( nw )
c
c local variables
      integer hgz, hgy, hgx, i, j, k, l, m, n, p, r, s, t
      real gm, rmax2, r2, rl2, x, y, z
c
      if( c3d )then
c set mass array to zero
        do j = 1, mesh( jgrid )
          grdmss( j, 1 ) = 0.
        end do
        rl2 = ( rtrunc( icmp ) * lscale )**2
c work over grid cells in 1st octant
        hgz = ngz / 2
        hgy = ngy / 2
        hgx = ngx / 2
        do k = 1, hgz
          z = k - 1
          p = ( k + hgz - 1 ) * ngxy
          do j = 1, hgy
            y = j - 1
            s = p + ( j + hgy - 1 ) * ngx
            do i = 1, hgx
              x = i - 1
              if( disc( icmp ) )then
                r2 = x * x + y * y
              else
                r2 = x * x + y * y + z * z
              end if
c confine to the largest inscribed sphere in the grid
              if( r2 .lt. rl2 )then
                call wcell( x, y, z, wt, rho )
c set pointers
                n = s + i + hgx
                l = n
                m = l + ngx
c distribute mass
                grdmss( l, 1 ) = grdmss( l, 1 ) + wt( 1 )
                grdmss( m, 1 ) = grdmss( m, 1 ) + wt( 3 )
                grdmss( l + 1, 1 ) = grdmss( l + 1, 1 ) + wt( 2 )
                grdmss( m + 1, 1 ) = grdmss( m + 1, 1 ) + wt( 4 )
                l = l + ngxy
                m = m + ngxy
                grdmss( l, 1 ) = grdmss( l, 1 ) + wt( 5 )
                grdmss( m, 1 ) = grdmss( m, 1 ) + wt( 7 )
                grdmss( l + 1, 1 ) = grdmss( l + 1, 1 ) + wt( 6 )
                grdmss( m + 1, 1 ) = grdmss( m + 1, 1 ) + wt( 8 )
              end if
            end do
          end do
        end do
c reflect about all 3 principal planes and add
c   n.b. wasteful use of a second array!
        n = 0
        l = -ngxy
        m = mesh( jgrid )
        do k = 1, ngz
          l = l + ngxy
          m = m - ngxy
          p = 0
          s = ngxy + 1
          do j = 1, ngy
            r = p + ngx + 1
            t = s - 1 - ngx
            do i = 1, ngx
              n = n + 1
              p = p + 1
              s = s - 1
              r = r - 1
              t = t + 1
              grdfld( n, 4 ) = grdmss( l + p, 1 ) + grdmss( m + p, 1 ) +
     +                         grdmss( l + s, 1 ) + grdmss( m + s, 1 ) +
     +                         grdmss( l + r, 1 ) + grdmss( m + r, 1 ) +
     +                         grdmss( l + t, 1 ) + grdmss( m + t, 1 )
            end do
          end do
        end do
c copy back to source array & rescale
        gm = lscale**3 * ts**2
        x = 0
        do i = 1, mesh( jgrid )
          x = x + grdfld( i, 4 )
          grdmss( i, 1 ) = gm * grdfld( i, 4 )
        end do
        print *, 'mass assigned to grid', x
c
      else if( p3d )then
c
        jrad( 1 ) = 1
        krad( 1 ) = nr( jgrid )
c note convention: x = r, y = theta, z = z
        do i = 1, mesh( jgrid )
          grdmss( i, 1 ) = 0.
        end do
        if( disc( icmp ) )then
          rmax2 = rgrid( jgrid )
        else
c restrict to largest inscribed sphere
          rmax2 = min( rgrid( jgrid ), zm( jgrid ) )**2
        end if
c work over grid radii
        y = 0.
        do i = 1, nr( jgrid ) - 1
          x = grofu( real( i - 1 ) )
c work over grid planes - assume reflection symmetry about z = 0
          do j = 1, ngz / 2
            z = dzg * real( j - 1 )
            r2 = x
            if( .not. disc( icmp ) )r2 = z * z + x * x
c points inside grid
            if( r2 .lt. rmax2 )then
              call wcell( x, y, z, wt, rho )
c set pointers for inner radius first
              l = ( ( i - 1 ) * ngz + ngz / 2 + j - 1 ) * na + 1
              m = l + na * ngz
c assign mass to one spoke only
              grdmss( l, 1 ) = grdmss( l, 1 ) + wt( 1 )
              grdmss( l + na, 1 ) = grdmss( l + na, 1 ) + wt( 5 )
              grdmss( m, 1 ) = grdmss( m, 1 ) + wt( 2 )
              grdmss( m + na, 1 ) = grdmss( m + na, 1 ) + wt( 6 )
            end if
          end do
        end do
c rotate in azimuth, reflect about z = 0 plane and rescale
        l = 0
        gm = lscale**3 * ts**2
        do i = 1, nr( jgrid )
          m = ( i * ngz - 1 ) * na
          do j = 1, ngz / 2 + 1
            y = gm * grdmss( m + 1, 1 )
            do k = 1, na
              grdmss( m + k, 1 ) = y
              l = l + 1
              grdmss( l, 1 ) = grdmss( l, 1 ) + y
            end do
            m = m - na
          end do
          l = l + na * ( ngz / 2 )
        end do
c
      else if( p3a )then
c
c note convention: x = r, z = z
        do i = 1, mesh( jgrid )
          grdmss( i, 1 ) = 0.
        end do
        if( disc( icmp ) )then
          rmax2 = rgrid( jgrid )
        else
c restrict to largest inscribed sphere
          rmax2 = min( rgrid( jgrid ), zm( jgrid ) )**2
        end if
        y = 0.
c work over grid planes - assume reflection symmetry about z = 0
        l = ( ngz / 2 ) * nr( jgrid )
        do j = 1, ngz / 2
          z = dzg * real( j - 1 )
c work over grid radii
          do i = 1, nr( jgrid )
            x = grofu( real( i - 1 ) )
            r2 = x
            if( .not. disc( icmp ) )r2 = z * z + x * x
            l = l + 1
c points inside grid
            if( r2 .lt. rmax2 )then
              call wcell( x, y, z, wt, rho )
              m = l + nr( jgrid )
c assign mass to one spoke only
              grdmss( l, 1 ) = grdmss( l, 1 ) + wt( 1 )
              grdmss( l + 1, 1 ) = grdmss( l + 1, 1 ) + wt( 2 )
              grdmss( m, 1 ) = grdmss( m, 1 ) + wt( 3 )
              grdmss( m + 1, 1 ) = grdmss( m + 1, 1 ) + wt( 4 )
            end if
          end do
        end do
c reflect about z = 0 plane and rescale
        l = 0
        m = mesh( jgrid ) - nr( jgrid )
        gm = lscale**3 * ts**2
        do j = 1, ngz / 2 + 1
          do i = 1, nr( jgrid )
            y = gm * grdmss( m + i, 1 )
            grdmss( m + i, 1 ) = y
            l = l + 1
            grdmss( l, 1 ) = grdmss( l, 1 ) + y
          end do
          m = m - nr( jgrid )
        end do
c
      else if( s3d )then
c
        if( s3lmax .gt. lmax )call crash( 'MASSDF',
     +                                    'wt array too small' )
c convert to grid units - int r**2 dr already includes lscale**3 factor
        gm = ts**2
c work over radii
        do i = 1, nr( jgrid ) - 1
          r2 = grofu( real( i - 1 ) )
c wcell returns the coeffs for the next radial cell
          call wcell( r2, 0., 0., wt, rho )
          j = 4 * s3ntm * ( i - 1 )
          do l = 0, s3lmax
            m = 2 * l
            j = j + 4 * l
c coefficients are purely real since m=0
            s3dmss( j + 1, 1, jgrid ) = gm * wt( m + 1 )
            s3dmss( j + 4 * s3ntm + 3, 1, jgrid ) = gm * wt( m + 2 )
          end do
        end do
      else
        call crash( 'MASSDF', 'Unrecognized grid' )
      end if
      return
      end
