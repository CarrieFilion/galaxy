      subroutine weight( jst, old )
c  Copyright (C) 2015, Jerry Sellwood
c
c Computes cell numbers, weights and radii for the current group of particles
c   from either the new or the old coordinates.  Particles off the grid are
c   flagged by setting the cell number to zero.
      implicit none
c
c calling arguments
      integer jst
      logical old
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c externals
      real grofu, uofgr
c
c local array
      real dp3( 3, 3 )
c
c local variables
      integer i, is, j, jt, ju, jx, jy, jz, k, l
      real dd, dt, du, dx, dy, dz, t, u, x, y, z
      include 'inc/pi.f'
c
c flag particles on the current grid
      do is = 1, jst
        nskip( is ) = ( iz( is ) .lt. nlists )
c set weights to zero in case they are not needed
        do ju = 1, 8
          wt( ju, is ) = 0
        end do
      end do
c add any particles on grid 2 to those to be skipped for acceleration
      if( old .and. ( jgrid .eq. 1 ) )then
        do is = 1, jst
          nskip( is ) = nskip( is ) .and. ( label( is ) .ne. 2 )
        end do
      end if
c add massless particles to those to be skipped for mass assignment
      if( .not. old )then
        do is = 1, jst
          i = iflag( is )
          nskip( is ) = nskip( is ) .and. ( .not. testp( i ) )
          nskip( is ) = nskip( is ) .and. ( label( is ) .eq. mhyb )
        end do
      end if
c
      if( p2d )then
        do is = 1, jst
          if( nskip( is ) )then
c get coordinates
            if( old )then
              x = oldc( 1, is ) - xcen( 1, jgrid )
              y = oldc( 2, is ) - xcen( 2, jgrid )
            else
              jz = max( iz( is ), 1 )
              x = newc( 1, is ) - xcpred( 1, jz, jgrid )
              y = newc( 2, is ) - xcpred( 2, jz, jgrid )
            end if
            rr( is ) = sqrt( x * x + y * y )
            u = uofgr( rr( is ) )
c            t = atan2( y, x ) / alpha
            t = atan2( dble( y ), dble( x ) ) *
     +                 dble( na ) / ( 2.d0 * pi )
            if( t .lt. 0. )t = t + real( na )
c ensure mass is assigned only to first sector
            if( ( .not. old ) .and.
     +          ( nsect .gt. 1 ) )t = mod( t, thmax )
c skip particles in central hole also
            if( ( rr( is ) .lt. rinh( jgrid ) ) )then
              nskip( is ) = .false.
            else
c compute cell number
              ju = u
              ju = min( ju, nr( jgrid ) - 2 )
              jt = t
              ncl( is ) = ju * na + jt + 1
c compute weights
              du = u - real( ju )
              dt = t - real( jt )
              dd = du * dt
              wt( 1, is ) = 1. - du - dt + dd
              wt( 2, is ) =           dt - dd
              wt( 3, is ) =      du      - dd
              wt( 4, is ) =                dd
            end if
          end if
        end do
      else if( c2d )then
        do is = 1, jst
c skip particles off grid
          if( nskip( is ) )then
c compute cell number in complete ( ngx x ngy ) grid plane
            if( old )then
              x = oldc( 1, is ) - xcen( 1, jgrid )
              y = oldc( 2, is ) - xcen( 2, jgrid )
            else
              jz = max( iz( is ), 1 )
              x = newc( 1, is ) - xcpred( 1, jz, jgrid )
              y = newc( 2, is ) - xcpred( 2, jz, jgrid )
            end if
            rr( is ) = sqrt( x * x + y * y )
            x = x + xm
            y = y + ym
            jy = y
            jx = x
            ncl( is ) = ( jy + 1 ) * ngx + jx + 2
c compute weights - CIC scheme
            dx = x - real( jx )
            dy = y - real( jy )
            wt( 1, is ) = ( 1. - dx ) * ( 1. - dy )
            wt( 2, is ) =        dx   * ( 1. - dy )
            wt( 3, is ) = ( 1. - dx ) *        dy
            wt( 4, is ) =        dx   *        dy
          end if
        end do
      else if( c3d )then
        do is = 1, jst
c skip particles off grid
          if( nskip( is ) )then
c compute cell number in complete ( ngx x ngy ) grid plane
            if( old )then
              x = oldc( 1, is ) + xm - xcen( 1, jgrid )
              y = oldc( 2, is ) + ym - xcen( 2, jgrid )
              z = oldc( 3, is ) + zm( jgrid ) - xcen( 3, jgrid )
            else
              x = newc( 1, is ) + xm - xcpred( 1, 1, jgrid )
              y = newc( 2, is ) + ym - xcpred( 2, 1, jgrid )
              z = newc( 3, is ) + zm( jgrid ) - xcpred( 3, 1, jgrid )
            end if
            jz = z
            jy = y
            jx = x
            ncl( is ) = ( jy + 1 ) * ngx + jx + 2
c compute weights - CIC scheme
            dx = x - real( jx )
            dy = y - real( jy )
            dz = z - real( jz )
            wt( 1, is ) = ( 1. - dx ) * ( 1. - dy ) * ( 1. - dz )
            wt( 2, is ) =        dx   * ( 1. - dy ) * ( 1. - dz )
            wt( 3, is ) = ( 1. - dx ) *        dy   * ( 1. - dz )
            wt( 4, is ) =        dx   *        dy   * ( 1. - dz )
            wt( 5, is ) = ( 1. - dx ) * ( 1. - dy ) *        dz
            wt( 6, is ) =        dx   * ( 1. - dy ) *        dz
            wt( 7, is ) = ( 1. - dx ) *        dy   *        dz
            wt( 8, is ) =        dx   *        dy   *        dz
          end if
        end do
      else if( p3a )then
        do is = 1, jst
c skip particles off grid
          if( nskip( is ) )then
c get coordinates
            if( old )then
              x = oldc( 1, is ) - xcen( 1, jgrid )
              y = oldc( 2, is ) - xcen( 2, jgrid )
              z = oldc( 3, is ) + zm( jgrid ) - xcen( 3, jgrid )
            else
              jz = max( iz( is ), 1 )
              x = newc( 1, is ) - xcpred( 1, jz, jgrid )
              y = newc( 2, is ) - xcpred( 2, jz, jgrid )
              z = newc( 3, is ) + zm( jgrid ) - xcpred( 3, jz, jgrid )
            end if
            u = sqrt( x * x + y * y )
            u = uofgr( u )
            z = z / dzg
c compute cell number
            ju = u
            jz = z
            jz = min( jz, ngz - 2 )
            ncl( is ) = jz * nr( jgrid ) + ju + 1
c compute weights
            du = u - real( ju )
            dz = z - real( jz )
            dd = du * dz
            wt( 1, is ) = 1. - du - dz + dd
            wt( 2, is ) =      du      - dd
            wt( 3, is ) =           dz - dd
            wt( 4, is ) =                dd
          end if
        end do
      else if( p3d )then
        if( jmass .eq. 2 )then
          do is = 1, jst
            if( nskip( is ) )then
c get coordinates
              if( old )then
                x = oldc( 1, is ) - xcen( 1, jgrid )
                y = oldc( 2, is ) - xcen( 2, jgrid )
                z = oldc( 3, is ) + zm( jgrid ) - xcen( 3, jgrid )
              else
                jz = max( iz( is ), 1 )
                x = newc( 1, is ) - xcpred( 1, jz, jgrid )
                y = newc( 2, is ) - xcpred( 2, jz, jgrid )
                z = newc( 3, is ) + zm( jgrid ) - xcpred( 3, jz, jgrid )
              end if
              rr( is ) = sqrt( x * x + y * y )
              u = uofgr( rr( is ) )
              if( rr( is ) .gt. 0. )then
c                t = atan2( y, x ) / alpha
                t = atan2( dble( y ), dble( x ) ) *
     +                     dble( na ) / ( 2.d0 * pi )
                if( t .lt. 0. )t = t + real( na )
              else
                t = 0
              end if
c ensure mass is assigned only to first sector
              if( ( .not. old ) .and.
     +            ( nsect .gt. 1 ) )t = mod( t, thmax )
c skip particles in central hole
              if( rr( is ) .lt. rinh( jgrid ) )then
                nskip( is ) = .false.
              else
c compute cell number
                ju = u
                ju = min( ju, nr( jgrid ) - 2 )
                jt = t
                z = z / dzg
                jz = z
                jz = max( jz, 0 )
                jz = min( jz, ngz - 2 )
                ncl( is ) = ( ju * ngz + jz ) * na + jt + 1
c compute weights
                du = u - real( ju )
                dt = t - real( jt )
                dz = z - real( jz )
                wt( 1, is ) = ( 1. - dt ) * ( 1. - dz ) * ( 1. - du )
                wt( 2, is ) =        dt   * ( 1. - dz ) * ( 1. - du )
                wt( 3, is ) = ( 1. - dt ) *        dz   * ( 1. - du )
                wt( 4, is ) =        dt   *        dz   * ( 1. - du )
                wt( 5, is ) = ( 1. - dt ) * ( 1. - dz ) *        du
                wt( 6, is ) =        dt   * ( 1. - dz ) *        du
                wt( 7, is ) = ( 1. - dt ) *        dz   *        du
                wt( 8, is ) =        dt   *        dz   *        du
              end if
            end if
          end do
        else if( jmass .eq. 3 )then
          do is = 1, jst
            if( nskip( is ) )then
c get coordinates
              if( old )then
                x = oldc( 1, is ) - xcen( 1, jgrid )
                y = oldc( 2, is ) - xcen( 2, jgrid )
                z = oldc( 3, is ) + zm( jgrid ) - xcen( 3, jgrid )
              else
                jz = max( iz( is ), 1 )
                x = newc( 1, is ) - xcpred( 1, jz, jgrid )
                y = newc( 2, is ) - xcpred( 2, jz, jgrid )
                z = newc( 3, is ) + zm( jgrid ) - xcpred( 3, jz, jgrid )
              end if
              rr( is ) = sqrt( x * x + y * y )
              u = uofgr( rr( is ) )
              if( rr( is ) .gt. 0. )then
c                t = atan2( y, x ) / alpha
                t = atan2( dble( y ), dble( x ) ) *
     +                     dble( na ) / ( 2.d0 * pi )
                if( t .lt. 0. )t = t + real( na )
              else
                t = 0
              end if
c ensure mass is assigned only to first sector
              if( ( .not. old ) .and.
     +            ( nsect .gt. 1 ) )t = mod( t, thmax )
c skip particles in central hole
              if( rr( is ) .lt. rinh( jgrid ) )then
                nskip( is ) = .false.
              else
c compute cell number - central cell of the 27
                ju = u + .5
                ju = min( ju, nr( jgrid ) - 2 )
                ju = max( ju, 1 )
                jt = t + .5
                z = z / dzg
                jz = z + .5
                jz = max( jz, 0 )
                jz = min( jz, ngz - 2 )
                ncl( is ) = ( ju * ngz + jz ) * na + jt + 1
                if( jt .eq. na )ncl( is ) = ncl( is ) - na
c compute weights
                do i = 1, 3
                  if( i .eq. 1 )dd = t + .5 - real( jt )
                  if( i .eq. 2 )dd = z + .5 - real( jz )
                  if( i .eq. 3 )dd = u + .5 - real( ju )
                  dp3( 1, i ) = 1. - 2. * dd + dd**2
                  dp3( 2, i ) = 1. + 2. * dd - 2. * dd**2
                  dp3( 3, i ) = dd**2
                end do
                l = 0
                do k = 1, 3
                  do j = 1, 3
                    do i = 1, 3
                      l = l + 1
                      wt3( l, is ) = .125 * dp3( i, 1 ) * dp3( j, 2 ) *
     +                                      dp3( k, 3 )
                    end do
                  end do
                end do
              end if
            end if
          end do
        end if
      else if( s3d )then
        do is = 1, jst
          if( old )nskip( is ) = nskip( is ) .or. offanal
c skip particles off grid
          if( nskip( is ) )then
c get coordinates
            if( old )then
              x = oldc( 1, is ) - xcen( 1, jgrid )
              y = oldc( 2, is ) - xcen( 2, jgrid )
              z = oldc( 3, is ) - xcen( 3, jgrid )
            else
              jz = max( iz( is ), 1 )
              x = newc( 1, is ) - xcpred( 1, jz, jgrid )
              y = newc( 2, is ) - xcpred( 2, jz, jgrid )
              z = newc( 3, is ) - xcpred( 3, jz, jgrid )
            end if
            rr( is ) = sqrt( x * x + y * y + z * z )
            if( old )then
c cell number and weights for accelerations
              u = uofgr( rr( is ) )
              ju = u + 1
              ncl( is ) = ju
              if( ju .ge. nr( jgrid ) )then
c particle outside outer shell
                ncl( is ) = nr( jgrid ) - 1
                wt( 1, is ) = 0
                wt( 2, is ) = 1
              else
c compute weights
                du = u - real( ju - 1 )
                wt( 1, is ) = 1. - du
                wt( 2, is ) =      du
              end if
            else
c cell number and weights for mass assignment - no interpolation
c              ju = uofgr( rr( is ) ) + 1
c              ncl( is ) = ju
c cell number and weights for mass assignment - linear interpolation
              u = uofgr( rr( is ) )
              ju = u + 1.5
              ncl( is ) = min( ju, nr( jgrid ) )
              if( ju .eq. 1 )then
c particle less than half a cell from the centre
c                wt( 1, is ) = 1. - u
c                wt( 2, is ) =      u
                wt( 1, is ) = 0
                wt( 2, is ) = 1
c              else if( ju .ge. nr( jgrid ) )then
c particle in outer shell
c                ncl( is ) = nr( jgrid )
c                wt( 1, is ) = 1
c                wt( 2, is ) = 0
              else
                du = u - real( ju - 1 )
                wt( 1, is ) = .5 - du
                wt( 2, is ) = .5 + du
c radii of mid-points of fragments
                wt( 3, is ) = grofu( .5 * ( real( ju - 1 ) + u - .5 ) )
                wt( 4, is ) = grofu( .5 * ( real( ju - 1 ) + u + .5 ) )
              end if
            end if
          end if
        end do
c SFP methods
      else if( sf2d .or. sf3d )then
c select limiting radius
        do is = 1, jst
          if( nskip( is ) )then
c jz for xcpred - not properly defined when recentering was introduced
            jz = max( iz( is ), 1 )
c find and store cylindrical radius
            if( old )then
              x = oldc( 1, is ) - xcen( 1, jgrid )
              y = oldc( 2, is ) - xcen( 2, jgrid )
            else
              x = newc( 1, is ) - xcpred( 1, jz, jgrid )
              y = newc( 2, is ) - xcpred( 2, jz, jgrid )
            end if
            rr( is ) = sqrt( x * x + y * y )
c plane number and weights
            if( sf3d )then
              if( old )then
                z =
     +    ( oldc( 3, is ) - xcpred( 3, jz, jgrid ) + zm( jgrid ) ) / dzg
              else
                z =
     +    ( newc( 3, is ) - xcpred( 3, jz, jgrid ) + zm( jgrid ) ) / dzg
              end if
c this jz is plane number
              jz = z
c particles outside range of grid planes - possible only for accelerations
              if( z .lt. 0. )then
                ncl( is ) = -1
                wt( 1, is ) =
     +                  -oldc( 3, is ) - zm( jgrid ) + xcen( 3, jgrid )
                wt( 2, is ) = 1
              else if( jz .gt. ngz - 2 )then
                ncl( is ) = ngz
                wt( 1, is ) =
     +                   oldc( 3, is ) - zm( jgrid ) + xcen( 3, jgrid )
                wt( 2, is ) = -1
              else
c particles within range of grid planes
                z = z - real( jz )
                ncl( is ) = jz + 1
                wt( 1, is ) = 1. - z
                wt( 2, is ) = z
              end if
            else
              ncl( is ) = 1
            end if
c flag particles outside active region
            if( ( rr( is ) .gt. rgrid( jgrid ) ) .or.
     +          ( hole .and. ( rr( is ) .lt. rinh( jgrid ) ) ) )then
              nskip( is ) = .false.
              rr( is ) = min( rr( is ), rgrid( jgrid ) )
            end if
          end if
        end do
        if( uqmass .and. ( .not. old ) )call crash( 'WEIGHT',
     +             'Individual particle masses not programmed for SFP' )
      else
        call crash( 'WEIGHT', 'Unrecognised method' )
      end if
c allow for unequal mass particles, if needed
      if( uqmass .and. ( .not. old ) )then
        jt = 8
        if( twod .or. p3a )jt = 4
        if( s3d )jt = 2
        do is = 1, jst
          do ju = 1, jt
            wt( ju, is ) = wt( ju, is ) * pwt( is )
          end do
        end do
      end if
      return
      end
