      subroutine difpot
c  Copyright (C) 2015, Jerry Sellwood
c
c written by Jerry Sellwood for the galaxy simulation code (c3d only)
c
c tabulates the three Cartesian acceleration components on two adjacent planes
c   of grid points
c
c The acceleration components are determined by differencing the pre-calculated
c   potential array.  In order to avoid storing each of these components at
c   every grid point, which would require three full grids of extra memory, the
c   potential differences are determined for just two planes at a time.  (Two
c   planes are needed for linear interpolation.)  This strategy requires all
c   particles in one grid plane to be accelerated before any others are
c   processed, and it is this requirement that forces the linked list
c   structure.
c
c The 10 point difference operator, devised by Andrew May and Richard James,
c   reduces force anisotropies well below those yielded by the obvious 2-point
c   formula.  [See Appendix A of Sellwood & Merritt (1994, ApJ v425, p530) for
c   a study of force quality.]
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variables
      integer i, j, l, m, n
      real f
c
      if( ( jplane .lt. 0 ) .or. ( jplane .gt. nplanes ) )then
        print *, 'jplane =', jplane
        call crash( 'DIFPOT', 'Nonsense jplane' )
      end if
      l = jplane * ngxy + ngx - 1
      n = ngx - 1
c work over plane
      do j = 2, ngy - 1
        l = l + 2
        n = n + 2
        do i = 2, ngx - 1
          l = l + 1
          m = l + ngxy
          n = n + 1
c x force components
          f = grdfld( l  + ngx + 1, 4 ) - grdfld( l  + ngx - 1, 4 ) +
     +        grdfld( l  - ngx + 1, 4 ) - grdfld( l  - ngx - 1, 4 ) +
     +        grdfld( l + ngxy + 1, 4 ) - grdfld( l + ngxy - 1, 4 ) +
     +        grdfld( l - ngxy + 1, 4 ) - grdfld( l - ngxy - 1, 4 ) +
     + 2. * ( grdfld( l        + 1, 4 ) - grdfld( l        - 1, 4 ) )
          c3dfld( n, 1 ) = -f / ( 12. * dh( 3 ) )
          f = grdfld( m  + ngx + 1, 4 ) - grdfld( m  + ngx - 1, 4 ) +
     +        grdfld( m  - ngx + 1, 4 ) - grdfld( m  - ngx - 1, 4 ) +
     +        grdfld( m + ngxy + 1, 4 ) - grdfld( m + ngxy - 1, 4 ) +
     +        grdfld( m - ngxy + 1, 4 ) - grdfld( m - ngxy - 1, 4 ) +
     + 2. * ( grdfld( m        + 1, 4 ) - grdfld( m        - 1, 4 ) )
          c3dfld( n + ngxy, 1 ) = -f / ( 12. * dh( 3 ) )
c y force components
          f =
     +       grdfld( l    + 1 + ngx, 4 ) - grdfld( l    + 1 - ngx, 4 ) +
     +       grdfld( l    - 1 + ngx, 4 ) - grdfld( l    - 1 - ngx, 4 ) +
     +       grdfld( l + ngxy + ngx, 4 ) - grdfld( l + ngxy - ngx, 4 ) +
     +       grdfld( l - ngxy + ngx, 4 ) - grdfld( l - ngxy - ngx, 4 ) +
     +2. * ( grdfld( l        + ngx, 4 ) - grdfld( l        - ngx, 4 ) )
          c3dfld( n, 2 ) = -f / ( 12. * dh( 2 ) )
          f =
     +       grdfld( m    + 1 + ngx, 4 ) - grdfld( m    + 1 - ngx, 4 ) +
     +       grdfld( m    - 1 + ngx, 4 ) - grdfld( m    - 1 - ngx, 4 ) +
     +       grdfld( m + ngxy + ngx, 4 ) - grdfld( m + ngxy - ngx, 4 ) +
     +       grdfld( m - ngxy + ngx, 4 ) - grdfld( m - ngxy - ngx, 4 ) +
     +2. * ( grdfld( m        + ngx, 4 ) - grdfld( m        - ngx, 4 ) )
          c3dfld( n + ngxy, 2 ) = -f / ( 12. * dh( 2 ) )
c z force components
          f =
     +       grdfld( l   + 1 + ngxy, 4 ) - grdfld( l   + 1 - ngxy, 4 ) +
     +       grdfld( l   - 1 + ngxy, 4 ) - grdfld( l   - 1 - ngxy, 4 ) +
     +       grdfld( l + ngx + ngxy, 4 ) - grdfld( l + ngx - ngxy, 4 ) +
     +       grdfld( l - ngx + ngxy, 4 ) - grdfld( l - ngx - ngxy, 4 ) +
     +2. * ( grdfld( l       + ngxy, 4 ) - grdfld( l       - ngxy, 4 ) )
          c3dfld( n, 3 ) = -f / ( 12. * dh( 1 ) )
          f =
     +       grdfld( m   + 1 + ngxy, 4 ) - grdfld( m   + 1 - ngxy, 4 ) +
     +       grdfld( m   - 1 + ngxy, 4 ) - grdfld( m   - 1 - ngxy, 4 ) +
     +       grdfld( m + ngx + ngxy, 4 ) - grdfld( m + ngx - ngxy, 4 ) +
     +       grdfld( m - ngx + ngxy, 4 ) - grdfld( m - ngx - ngxy, 4 ) +
     +2. * ( grdfld( m       + ngxy, 4 ) - grdfld( m       - ngxy, 4 ) )
          c3dfld( n + ngxy, 3 ) = -f / ( 12. * dh( 1 ) )
        end do
      end do
      return
      end
