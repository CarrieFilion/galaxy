      subroutine grdacc( jst )
c  Copyright (C) 2015, Jerry Sellwood
c
c Returns acceleration acting on each particle in the current group arising
c   from contributions on the grid only.  The acceleration components are the
c   weighted sums of values the nearest cell corners.
c The appropriate weights are obtained by a call to subroutine WEIGHT, which
c   assumes that the "old" coordinates contain the appropriate position.
c   Linear interpolation between mesh points is hard-wired except for p3d.
c No supplementary or external components are added in this routine.
      use aarrays
      implicit none
c
c calling argument
      integer jst
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
c external
      real uofgr
c
c local arrays
      integer ic( 27 )
      real ac( 2 )
c
c local variables
      integer i, is, j
      logical end
      real rad, u, x, y
c
c compute grid cell number and weights for interpolation
      call weight( jst, .true. )
c
      if( p2d )then
        do is = 1, jst
          if( nskip( is ) )then
            u = uofgr( rr( is ) )
c tabulate locations of cell corners in each array
            end = mod( ncl( is ), na ) .eq. 0
            ic( 1 ) = ncl( is )
            ic( 3 ) = ncl( is ) + na
            ic( 2 ) = ncl( is ) + 1
            if( end )ic( 2 ) = ic( 2 ) - na
            ic( 4 ) = ic( 2 ) + na
c evaluate accelerations
            do i = 1, 2
              do j = 1, 4
                acc( i, is ) =
     +                 acc( i, is ) + grdfld( ic( j ), i ) * wt( j, is )
              end do
            end do
c potential value if needed
            if( phys )then
              do i = 1, 4
                gpot( is ) =
     +                   gpot( is ) + grdfld( ic( i ), 3 ) * wt( i, is )
              end do
c subtract self-energy
c              x = oldc( 1, is ) - xcen( 1, 1 )
c              y = oldc( 2, is ) - xcen( 2, 1 )
c              u = sqrt( x * x + y * y )
c              u = uofgr( u )
c              i = int( u ) + 1
c              u = wt( 1, is ) * wt( 2, is ) * self( 2, i ) +
c     +            wt( 3, is ) * wt( 4, is ) * self( 2, i + 1 )
c              u = u + ( wt( 1, is ) * wt( 3, is ) +
c     +                  wt( 2, is ) * wt( 4, is ) ) * self( 3, i ) +
c     +                ( wt( 1, is ) * wt( 4, is ) +
c     +                  wt( 2, is ) * wt( 3, is ) ) * self( 4, i )
c              u = 2. * u +
c     +            ( wt( 1, is )**2 + wt( 2, is )**2 ) * self( 1, i ) +
c     +            ( wt( 3, is )**2 + wt( 4, is )**2 ) * self( 1, i + 1 )
c              gpot( is ) = gpot( is ) - u
            end if
          end if
        end do
c
      else if( c2d )then
        do is = 1, jst
          if( nskip( is ) )then
c tabulate locations of cell corners in each array
            ic( 1 ) = ncl( is )
            ic( 2 ) = ncl( is ) + 1
            ic( 3 ) = ncl( is ) + ngx
            ic( 4 ) = ncl( is ) + ngx + 1
c evaluate accelerations
            do i = 1, 2
              do j = 1, 4
                acc( i, is ) =
     +                 acc( i, is ) + grdfld( ic( j ), i ) * wt( j, is )
              end do
            end do
c potential value if needed
            if( phys )then
              do i = 1, 4
                gpot( is ) =
     +                   gpot( is ) + grdfld( ic( i ), 3 ) * wt( i, is )
              end do
c subtract self-energy
c              x = dh( 1 )
c              y = dh( 2 )
c              u =
c     +         ( wt( 1, is ) * wt( 3, is ) + wt( 2, is ) * wt( 4, is ) )
c     +                                           / sqrt( softl2 + y**2 )
c     +       + ( wt( 1, is ) * wt( 2, is ) + wt( 3, is ) * wt( 4, is ) )
c     +                                           / sqrt( softl2 + x**2 )
c     +       + ( wt( 2, is ) * wt( 3, is ) + wt( 1, is ) * wt( 4, is ) )
c     +                                    / sqrt( softl2 + x**2 + y**2 )
c              u = 2. * u + ( wt( 1, is )**2 + wt( 2, is )**2 +
c     +                       wt( 3, is )**2 + wt( 4, is )**2 ) / softl
c              u = u * pmass
c              gpot( is ) = gpot( is ) + u
            end if
          end if
        end do
c
      else if( c3d )then
        do is = 1, jst
          if( nskip( is ) )then
            if( ipln( is ) .ne. jplane )then
              jplane = ipln( is )
              call difpot
            end if
c tabulate locations of cell corners in each plane
            ic( 1 ) = ncl( is )
            ic( 2 ) = ncl( is ) + 1
            ic( 3 ) = ncl( is ) + ngx
            ic( 4 ) = ncl( is ) + ngx + 1
            do i = 1, 4
              ic( i + 4 ) = ic( i ) + ngxy
            end do
c evaluate accelerations
            do i = 1, 3
              do j = 1, 8
                acc( i, is ) =
     +                 acc( i, is ) + c3dfld( ic( j ), i ) * wt( j, is )
              end do
            end do
c potential value if needed
            if( phys )then
              do i = 1, 8
                j = ic( i ) + jplane * ngxy
                gpot( is ) = gpot( is ) + grdfld( j, 4 ) * wt( i, is )
              end do
            end if
          end if
        end do
c
      else if( p3a )then
        do is = 1, jst
          if( nskip( is ) )then
c tabulate locations of cell corners in each array
            ic( 1 ) = ncl( is )
            ic( 2 ) = ncl( is ) + 1
            ic( 3 ) = ncl( is ) + nr( jgrid )
            ic( 4 ) = ncl( is ) + nr( jgrid ) + 1
c evaluate accelerations
            do i = 1, 2
              ac( i ) = 0
              do j = 1, 4
                ac( i ) = ac( i ) + grdfld( ic( j ), i ) * wt( j, is )
              end do
            end do
c resolve radial acceleration into Cartesian components
            x = oldc( 1, is ) - xcen( 1, 1 )
            y = oldc( 2, is ) - xcen( 2, 1 )
            rad = sqrt( x * x + y * y )
            acc( 1, is ) = acc( 1, is ) + ac( 1 ) * x / rad
            acc( 2, is ) = acc( 2, is ) + ac( 1 ) * y / rad
            acc( 3, is ) = acc( 3, is ) + ac( 2 )
c potential value if needed
            if( phys )then
              do i = 1, 4
                gpot( is ) =
     +                   gpot( is ) + grdfld( ic( i ), 3 ) * wt( i, is )
              end do
            end if
          end if
        end do
c
      else if( p3d )then
        if( jmass .eq. 2 )then
          do is = 1, jst
            if( nskip( is ) )then
c tabulate locations of cell corners in each array
              end = mod( ncl( is ), na ) .eq. 0
              ic( 1 ) = ncl( is )
              ic( 3 ) = ncl( is ) + na
              ic( 2 ) = ncl( is ) + 1
              if( end )ic( 2 ) = ic( 2 ) - na
              ic( 4 ) = ic( 2 ) + na
              ic( 5 ) = ncl( is ) + ngz * na
              ic( 7 ) = ic( 5 ) + na
              ic( 6 ) = ic( 5 ) + 1
              if( end )ic( 6 ) = ic( 6 ) - na
              ic( 8 ) = ic( 6 ) + na
c evaluate accelerations
              do i = 1, 3
                do j = 1, 8
                  acc( i, is ) =
     +                 acc( i, is ) + grdfld( ic( j ), i ) * wt( j, is )
                end do
              end do
c potential value if needed
              if( phys )then
                do i = 1, 8
                  gpot( is ) =
     +                   gpot( is ) + grdfld( ic( i ), 4 ) * wt( i, is )
                end do
c subtract self-energy
c                x = oldc( 1, is ) - xcen( 1, 1 )
c                y = oldc( 2, is ) - xcen( 2, 1 )
c                u = sqrt( x * x + y * y )
c                u = uofgr( u )
c                i = int( u ) + 1
c                u = ( wt( 1, is ) * wt( 2, is ) +
c     +                wt( 3, is ) * wt( 4, is ) ) * self( 2, i ) +
c     +              ( wt( 5, is ) * wt( 6, is ) +
c     +                wt( 7, is ) * wt( 8, is ) ) * self( 2, i + 1 )
c                u = u +
c     +              ( wt( 1, is ) * wt( 5, is ) +
c     +                wt( 2, is ) * wt( 6, is ) +
c     +                wt( 3, is ) * wt( 7, is ) +
c     +                wt( 4, is ) * wt( 8, is ) ) * self( 3, i )
c                u = u +
c     +              ( wt( 1, is ) * wt( 6, is ) +
c     +                wt( 2, is ) * wt( 5, is ) +
c     +                wt( 3, is ) * wt( 8, is ) +
c     +                wt( 4, is ) * wt( 7, is ) ) * self( 4, i )
c                u = u +
c     +              ( wt( 1, is ) * wt( 3, is ) +
c     +                wt( 2, is ) * wt( 4, is ) ) * self( 7, i ) +
c     +              ( wt( 5, is ) * wt( 7, is ) +
c     +                wt( 6, is ) * wt( 8, is ) ) * self( 7, i + 1 )
c                u = u +
c     +              ( wt( 1, is ) * wt( 4, is ) +
c     +                wt( 2, is ) * wt( 3, is ) ) * self( 8, i ) +
c     +              ( wt( 5, is ) * wt( 8, is ) +
c     +                wt( 6, is ) * wt( 7, is ) ) * self( 8, i + 1 )
c                u = u +
c     +              ( wt( 1, is ) * wt( 7, is ) +
c     +                wt( 2, is ) * wt( 8, is ) +
c     +                wt( 3, is ) * wt( 5, is ) +
c     +                wt( 4, is ) * wt( 6, is ) ) * self( 9, i )
c                u = u +
c     +              ( wt( 1, is ) * wt( 8, is ) +
c     +                wt( 2, is ) * wt( 7, is ) +
c     +                wt( 3, is ) * wt( 6, is ) +
c     +               wt( 4, is ) * wt( 5, is ) ) * self( 10, i )
c                u = 2 * u +
c     +            ( wt( 1, is )**2 + wt( 2, is )**2 +
c     +              wt( 5, is )**2 + wt( 6, is )**2 ) * self( 1, i ) +
c     +            ( wt( 3, is )**2 + wt( 4, is )**2 +
c     +              wt( 7, is )**2 + wt( 8, is )**2 ) * self( 1, i + 1 )
c                gpot( is ) = gpot( is ) - u
              end if
            end if
          end do
        else if( jmass .eq. 3 )then
          do is = 1, jst
            if( nskip( is ) )then
c tabulate locations of cell corners in each array
              ic( 14 ) = ncl( is )
              ic( 15 ) = ncl( is ) + 1
              ic( 13 ) = ncl( is ) - 1
              if( mod( ncl( is ), na ) .eq. 0 )ic( 15 ) = ic( 15 ) - na
              if( mod( ncl( is ), na ) .eq. 1 )ic( 13 ) = ic( 13 ) + na
              do j = 13, 15
                ic( j - 3 ) = ic( j ) - na
                ic( j + 3 ) = ic( j ) + na
              end do
              do j = 10, 18
                ic( j - 9 ) = ic( j ) - ngz * na
                ic( j + 9 ) = ic( j ) + ngz * na
              end do
c evaluate accelerations
              do i = 1, 3
                do j = 1, 27
                  acc( i, is ) =
     +                acc( i, is ) + grdfld( ic( j ), i ) * wt3( j, is )
                end do
              end do
c potential value if needed
              if( phys )then
                do i = 1, 27
                  gpot( is ) =
     +                  gpot( is ) + grdfld( ic( i ), 4 ) * wt3( i, is )
                end do
              end if
            end if
          end do
        end if
c
      else
        call crash( 'GRDACC', 'Unrecognised grid' )
      end if
      return
      end
