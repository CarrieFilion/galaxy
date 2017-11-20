      subroutine massgn( jst )
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to assign the masses of the current group of particles to the grid
c   or to increment the expansion coefficients for SFP methods.
c
c It is necessary to test for particles off the grid, since this routine
c   is normally called after the particles have been stepped forward, and
c   they could have left the grid
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
      include 'inc/bdpps.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c local array
      integer ic( 27 )
c
c local variables
      integer is, j, jz, k, l, n
      logical end
c
c      if( ( jgrid .gt. 1 ) .and. ( .not. s3d ) )call crash( 'MASSGN',
c     +                                         'Option not programmed' )
c SCF method
      if( scf )then
        call scfadd( jst )
      else if( dr3d )then
c direct-N method
        do is = 1, jst
          if( label( is ) .eq. jgrid )then
            n = 1 + loc( is ) / nwpp - kdrct
            do j = 1, 3
              drpt( j + 4, n ) = newc( j, is )
            end do
          end if
        end do
      else if( .not. bht )then
c determine cell numbers and weights
        call weight( jst, .false. )
c SFP methods
        if( sf2d .or. sf3d )then
c update newmaxr
          do is = 1, jst
            newmaxr = max( rr( is ), newmaxr )
          end do
          if( ( basset .eq. 'bess' ) .or. ( basset .eq. 'lgsp' ) )then
            call sfpadd( jst )
          else if( basset .ne. 'ablj' )then
            call crash( 'MASSGN', 'Unrecognized basis' )
          end if
        else if( s3d )then
c PM+SH method
          call s3dsum( jst )
c grid methods
        else if( p2d )then
          do is = 1, jst
c skip particles off this grid
            if( nskip( is ) )then
c determine radial range of distributed mass
              j = ( ncl( is ) - 1 ) / na
              jz = max( iz( is ), 1 )
              jrad( jz ) = min( jrad( jz ), j + 1 )
              krad( jz ) = max( krad( jz ), j + 2 )
c set pointers
              j = ncl( is )
              k = ncl( is ) + na
c distribute mass
              grdmss( j, jz ) = grdmss( j, jz ) + wt( 1, is )
              grdmss( k, jz ) = grdmss( k, jz ) + wt( 3, is )
              j = j + 1
              if( mod( ncl( is ), na ) .eq. 0 )j = j - na
              k = j + na
              grdmss( j, jz ) = grdmss( j, jz ) + wt( 2, is )
              grdmss( k, jz ) = grdmss( k, jz ) + wt( 4, is )
            end if
          end do
        else if( c2d )then
          do is = 1, jst
c skip particles off this grid
            if( nskip( is ) )then
c set pointers
              jz = max( iz( is ), 1 )
              j = ncl( is )
              k = ncl( is ) + ngx
c distribute mass
              grdmss( j, jz ) = grdmss( j, jz ) + wt( 1, is )
              grdmss( k, jz ) = grdmss( k, jz ) + wt( 3, is )
              grdmss( j + 1, jz ) = grdmss( j + 1, jz ) + wt( 2, is )
              grdmss( k + 1, jz ) = grdmss( k + 1, jz ) + wt( 4, is )
            end if
          end do
        else if( c3d )then
          do is = 1, jst
c skip particles off this grid
            if( nskip( is ) )then
c set pointers
              j = ngxy * ipln( is ) + ncl( is )
              k = j + ngx
c distribute mass
              grdmss( j, 1 ) = grdmss( j, 1 ) + wt( 1, is )
              grdmss( k, 1 ) = grdmss( k, 1 ) + wt( 3, is )
              grdmss( j + 1, 1 ) = grdmss( j + 1, 1 ) + wt( 2, is )
              grdmss( k + 1, 1 ) = grdmss( k + 1, 1 ) + wt( 4, is )
              j = j + ngxy
              k = k + ngxy
              grdmss( j, 1 ) = grdmss( j, 1 ) + wt( 5, is )
              grdmss( k, 1 ) = grdmss( k, 1 ) + wt( 7, is )
              grdmss( j + 1, 1 ) = grdmss( j + 1, 1 ) + wt( 6, is )
              grdmss( k + 1, 1 ) = grdmss( k + 1, 1 ) + wt( 8, is )
            end if
          end do
        else if( p3a )then
          do is = 1, jst
c skip particles off this grid
            if( nskip( is ) )then
c set pointers
              jz = max( iz( is ), 1 )
              j = ncl( is )
              k = ncl( is ) + nr( jgrid )
c distribute mass
              grdmss( j, jz ) = grdmss( j, jz ) + wt( 1, is )
              grdmss( k, jz ) = grdmss( k, jz ) + wt( 3, is )
              grdmss( j + 1, jz ) = grdmss( j + 1, jz ) + wt( 2, is )
              grdmss( k + 1, jz ) = grdmss( k + 1, jz ) + wt( 4, is )
            end if
          end do
        else if( p3d )then
          if( jmass .eq. 2 )then
            do is = 1, jst
c skip particles off this grid
              if( nskip( is ) )then
c determine radial range of distributed mass
                j = ( ncl( is ) - 1 ) / ( na * ngz )
                jz = max( iz( is ), 1 )
                jrad( jz ) = min( jrad( jz ), j + 1 )
                krad( jz ) = max( krad( jz ), j + 2 )
c set pointers
                j = ncl( is )
                l = j + ngz * na
c distribute mass
                k = j + na
                grdmss( j, jz ) = grdmss( j, jz ) + wt( 1, is )
                grdmss( k, jz ) = grdmss( k, jz ) + wt( 3, is )
                j = j + 1
                end = mod( ncl( is ), na ) .eq. 0
                if( end )j = j - na
                k = j + na
                grdmss( j, jz ) = grdmss( j, jz ) + wt( 2, is )
                grdmss( k, jz ) = grdmss( k, jz ) + wt( 4, is )
                k = l + na
                grdmss( l, jz ) = grdmss( l, jz ) + wt( 5, is )
                grdmss( k, jz ) = grdmss( k, jz ) + wt( 7, is )
                l = l + 1
                if( end )l = l - na
                k = l + na
                grdmss( l, jz ) = grdmss( l, jz ) + wt( 6, is )
                grdmss( k, jz ) = grdmss( k, jz ) + wt( 8, is )
              end if
            end do
          else if( jmass .eq. 3 )then
            do is = 1, jst
c skip particles off this grid
              if( nskip( is ) )then
c determine radial range of distributed mass
                j = ( ncl( is ) - 1 ) / ( na * ngz )
                jz = max( iz( is ), 1 )
                jrad( jz ) = min( jrad( jz ), j  )
                krad( jz ) = max( krad( jz ), j + 2 )
c tabulate locations of cell corners
                ic( 14 ) = ncl( is )
                ic( 15 ) = ic( 14 ) + 1
                ic( 13 ) = ic( 14 ) - 1
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
c assign mass
                do j = 1, 27
                  k = ic( j )
                  grdmss( k, jz ) = grdmss( k, jz ) + wt3( j, is )
                end do
              end if
            end do
          end if
        else
          call crash( 'MASSGN', 'Unrecognised method' )
        end if
      end if
      return
      end
