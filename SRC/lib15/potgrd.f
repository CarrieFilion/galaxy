      real function potgrd( xf, yf, zf )
c  Copyright (C) 2015, Jerry Sellwood
c
c returns potential on the grid at an arbitrary field point by computing
c   the work required from the centre.  The work required is determined
c   by computing a line integral of the grid force starting from centre
c   and moving outwards in directions parallel to the mesh axes.
c
c   3-fold reflection symmetry assumed for c3d
c   Axial symmetry assumed for p3d
      use aarrays
      implicit none
c
c calling arguments - position of field point in mesh spaces
      real xf, yf, zf
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      real aold, potold, zold
      logical recalc
      common / potsav / recalc, aold, zold, potold
c
c externals
      real grofu, uofgr
      real*8 quad_tab
c
c local allocatable arrays
      real, allocatable :: f(:,:,:,:)
      real*8, allocatable :: abscis(:)
      real*8, allocatable :: weight(:)
      save abscis, f, weight
c
c local variables
      integer i, iuse, ix, iy, jz, j, k, n
      logical xfst, yfst, zfst
      real a, dx, dy, dz, fx1, fx2, fy1, fy2, fz1, fz2
      real w, x, x1, x2, y, y1, y2, z, z1, z2
      real*8 r1, r2
      save iuse
      data iuse / 0 /
c
c tabulate forces
      if( c3d )then
        if( recalc )then
          if( iuse .eq. 0 )then
            allocate ( f( 3, ngx / 2, ngy / 2, ngz / 2 ) )
            iuse = 1
          end if
          do jz = 1, nplanes / 2
            jplane = jz + nplanes / 2 - 1
            call difpot
            do iy = 1, ngy / 2
              n = ngx * ( ngy / 2 + iy - 1 ) + ngx / 2
              do ix = 1, ngx / 2
                n = n + 1
                do i = 1, 3
                  f( i, ix, iy, jz ) = c3dfld( n, i )
                end do
              end do
            end do
          end do
          recalc = .false.
        end if
c impose reflection symmetry
        x = abs( xf )
        y = abs( yf )
        z = abs( zf )
c
        if( ( x .gt. xm ) .or. ( y .gt. ym ) .or.
     +      ( z .gt. zm( jgrid ) ) )then
          print *, xf, yf, zf
          call crash( 'POTGRD', 'Point outside grid' )
        end if
c form line integral along mesh axes starting from the centre
c   where forces are zero
        potgrd = 0.
        a = max( x, y, z )
        if( a .gt. 0. )then
c step first along longest axis
          xfst = a .eq. x
          yfst = ( a .eq. y ) .and. ( .not. xfst )
          zfst = ( a .eq. z ) .and. ( .not. ( xfst .or. yfst ) )
c
          x2 = 0
          y2 = 0
          z2 = 0
          ix = 0
          iy = 0
          jz = 0
          fx2 = 0
          fy2 = 0
          fz2 = 0
c work over line segments until done
          do while ( x - x2 + y - y2 + z - z2 .gt. 1.e-5 )
c determine next step
            x1 = x2
            y1 = y2
            z1 = z2
    2       if( xfst )then
              x2 = ix + 1
              x2 = min( x2, x )
              if( x1 .eq. x )then
                xfst = .false.
                a = max( y - y1, z - z1 )
                yfst = ( a .eq. y ) .and. ( y1 .lt. y )
                if( .not. yfst )zfst = a .eq. z
              end if
            end if
    3       if( yfst )then
              y2 = iy + 1
              y2 = min( y2, y )
              if( y1 .eq. y )then
                yfst = .false.
                a = max( x - x1, z - z1 )
                xfst = ( a .eq. x ) .and. ( x1 .lt. x )
                if( xfst )go to 2
                zfst = a .eq. z
              end if
            end if
            if( zfst )then
              z2 = jz + 1
              z2 = min( z2, z )
              if( z1 .eq. z )then
                zfst = .false.
                a = max( x - x1, y - y1 )
                xfst = ( a .eq. x ) .and. ( x1 .lt. x )
                if( xfst )go to 2
                yfst = a .eq. y
                if( yfst )go to 3
              end if
            end if
c evaluate force components at each end of this line segment
            fx1 = fx2
            fy1 = fy2
            fz1 = fz2
            ix = x2
            dx = x2 - real( ix )
            iy = y2
            dy = y2 - real( iy )
            jz = z2
            dz = z2 - real( jz )
            fx2 = 0
            fy2 = 0
            fz2 = 0
            do k = jz + 1, jz + 2
              dz = 1. - dz
              do j = iy + 1, iy + 2
                dy = 1. - dy
                do i = ix + 1, ix + 2
                  dx = 1. - dx
                  w = dx * dy * dz
                  if( w .gt. 0. )then
                    fx2 = fx2 + w * f( 1, i, j, k )
                    fy2 = fy2 + w * f( 2, i, j, k )
                    fz2 = fz2 + w * f( 3, i, j, k )
                  end if
                end do
              end do
            end do
c compute work done in grid units
            potgrd = potgrd - .5 * ( fx1 + fx2 ) * ( x2 - x1 )
     +                      - .5 * ( fy1 + fy2 ) * ( y2 - y1 )
     +                      - .5 * ( fz1 + fz2 ) * ( z2 - z1 )
          end do
        end if
c
      else if( p3d )then
c
        if( ( xf .gt. rgrid( jgrid ) ) .or.
     +      ( zf .gt. zm( jgrid ) ) )then
          print *, xf, yf, zf
          call crash( 'POTGRD', 'Point outside grid' )
        end if
        potgrd = 0.
c radial direction first to nearest interior ring
        if( xf .gt. 0. )then
          ix = uofgr ( xf )
          j = na * ( ngz / 2 ) + 1
          fx2 = grdfld( j, 1 )
          x2 = grofu( 0. )
          do i = 1, ix
            j = j + na * ngz
            fx1 = fx2
            fx2 = grdfld( j, 1 )
            x1 = x2
            x2 = grofu( real( i ) )
            potgrd = potgrd - .5 * ( fx1 + fx2 ) * ( x2 - x1 )
          end do
c sub-grid radial part
          dx = xf - x2
          j = j + na * ngz
          fx1 = fx2
          fx2 = grdfld( j, 1 )
          x1 = x2
          x2 = grofu( real( ix + 1 ) )
          fx2 = fx1 + dx * ( fx2 - fx1 ) / ( x2 - x1 )
          potgrd = potgrd - .5 * dx * ( fx1 + fx2 )
          dx = dx / ( x2 - x1 )
        else
          ix = 0
          dx = 0
        end if
c vertical direction next, to plane just below field point
        if( zf .ne. 0. )then
          jz = abs( zf ) / dzg
          j = ( ix * ngz + ngz / 2 ) * na + 1
          fz2 = ( 1. - dx ) * grdfld( j, 3 ) +
     +                 dx   * grdfld( j + na * ngz, 3 )
          do k = 1, jz
            j = j + na
            fz1 = fz2
            fz2 = ( 1.- dx ) * grdfld( j, 3 ) +
     +                  dx   * grdfld( j + na * ngz, 3 )
            potgrd = potgrd - .5 * ( fz1 + fz2 ) * dzg
          end do
c sub-grid vertical part
          dz = abs( zf ) - real( jz ) * dzg
          j = j + na
          fz1 = fz2
          fz2 = ( 1. - dx ) * grdfld( j, 3 ) +
     +                 dx   * grdfld( j + na * ngz, 3 )
          fz2 = fz1 + dz * ( fz2 - fz1 ) / dzg
          potgrd = potgrd - .5 * dz * ( fz1 + fz2 )
        end if
c
      else if( p3a )then
c
        if( ( xf .gt. rgrid( jgrid ) ) .or.
     +      ( zf .gt. zm( jgrid ) ) )then
          print *, xf, yf, zf
          call crash( 'POTGRD', 'Point outside grid' )
        end if
        potgrd = 0.
c radial direction first to nearest interior ring
        if( xf .gt. 0. )then
          ix = uofgr( xf )
          j = nr( jgrid ) * ( ngz / 2 ) + 1
          fx2 = grdfld( j, 1 )
          x2 = grofu( 0. )
          do i = 1, ix
            j = j + 1
            fx1 = fx2
            fx2 = grdfld( j, 1 )
            x1 = x2
            x2 = grofu( real( i ) )
            potgrd = potgrd - .5 * ( fx1 + fx2 ) * ( x2 - x1 )
          end do
c sub-grid radial part
          dx = xf - x2
          j = j + 1
          fx1 = fx2
          fx2 = grdfld( j, 1 )
          x1 = x2
          x2 = grofu( real( ix + 1 ) )
          fx2 = fx1 + dx * ( fx2 - fx1 ) / ( x2 - x1 )
          potgrd = potgrd - .5 * dx * ( fx1 + fx2 )
          dx = dx / ( x2 - x1 )
        else
          ix = 0
          dx = 0
        end if
c vertical direction next, to plane just below field point
        if( zf .ne. 0. )then
          jz = abs( zf ) / dzg
          j = nr( jgrid ) * ( ngz / 2 ) + ix + 1
          fz2 = ( 1. - dx ) * grdfld( j, 3 ) +
     +                 dx   * grdfld( j + 1, 3 )
          do k = 1, jz
            j = j + nr( jgrid )
            fz1 = fz2
            fz2 = ( 1. - dx ) * grdfld( j, 3 ) +
     +                   dx   * grdfld( j + 1, 3 )
            potgrd = potgrd - .5 * ( fz1 + fz2 ) * dzg
          end do
c sub-grid vertical part
          dz = abs( zf ) - real( jz ) * dzg
          j = j + nr( jgrid )
          fz1 = fz2
          fz2 = ( 1. - dx ) * grdfld( j, 3 ) +
     +                 dx   * grdfld( j + 1, 3 )
          fz2 = fz1 + dz * ( fz2 - fz1 ) / dzg
          potgrd = potgrd - .5 * dz * ( fz1 + fz2 )
        end if
c
      else
c
c generic method
        izone = 1
        potgrd = 0
        if( .not. allocated( abscis ) )then
          allocate( abscis( mbuff ) )
          allocate( weight( mbuff ) )
        end if
        a = sqrt( xf**2 + yf**2 )
        if( a .ne. aold )then
c radial part
          if( a .gt. 0. )then
            dx = rgrid( ngrid ) / real( mbuff - 1 )
c divide distance into equal segements
            n = a / dx + 1
            n = max( n, 2 )
            dx = a / real( n - 1 )
            x = a
            do i = 1, n - 1
c create fictitious particles
              abscis( i ) = x
              oldc( 1, i ) = x * xf / a
              oldc( 2, i ) = x * yf / a
              oldc( 3, i ) = 0
              iz( i ) = 1
              label( i ) = ngrid
              iflag( i ) = 2
              x = x - dx
            end do
c get the accelerations for these particles
            call getacc( n - 1 )
c sum two force components along polar radius vector
            do i = 1, n - 1
              weight( i ) = ( oldc( 1, i ) * acc( 1, i ) +
     +                        oldc( 2, i ) * acc( 2, i ) ) / abscis( i )
            end do
c point at the centre
            abscis( n ) = 0
            weight( n ) = 0
c compute work done
            if( n .gt. 3 )then
              r1 = quad_tab( abscis, weight, n, r2 )
              potgrd = r1 - r2
            else
c trapezium rule will do
              do i = 2, n
                potgrd = potgrd +
     +                    .5 * ( weight( i ) + weight( i - 1 ) ) *
     +                         ( abscis( i ) - abscis( i - 1 ) )
              end do
            end if
          end if
          aold = a
          potold = potgrd
        else
          potgrd = potold
        end if
c vertical part
        if( zf .gt. 0. )then
          dz = zm( ngrid ) / real( mbuff - 1 )
c divide distance into equal segements
          n = zf / dz + 1
          n = max( 2, n )
          dz = zf / real( n - 1 )
          z = zf
c divide distance into equal segements
          do i = 1, n - 1
c create fictitious particles
            abscis( i ) = z
            oldc( 1, i ) = xf
            oldc( 2, i ) = yf
            oldc( 3, i ) = z
            iz( i ) = 1
            label( i ) = ngrid
            iflag( i ) = 2
            z = z - dz
          end do
c get the accelerations for these particles
          call getacc( n - 1 )
          do i = 1, n - 1
            weight( i ) = acc( 3, i )
          end do
c point in the mid-plane - assume vertical force is zero there
          abscis( n ) = 0
          weight( n ) = 0
c compute work done
          if( n .gt. 3 )then
            r1 = quad_tab( abscis, weight, n, r2 )
            potgrd = potgrd + r1 - r2
          else
c trapezium rule will do
            do i = 2, n
              potgrd = potgrd +
     +                   .5 * ( weight( i ) + weight( i - 1 ) ) *
     +                        ( abscis( i ) - abscis( i - 1 ) )
            end do
          end if
        end if
      end if
      return
      end
