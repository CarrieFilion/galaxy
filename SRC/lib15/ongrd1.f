      logical function ongrd1( is )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Determines whether particle is inside the smaller grid of a hybrid code
c   on the basis of its newc position for APLF time stepping or the
c   oldc position otherwise.  This function is irrelevant for a simple
c   (non-hybrid) grid, as it is just the negation of OFFGRD
c
c Called from REZGRP
c
c calling argument
      integer is
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
c "local" common
      integer iuse
      logical lc2d, lc3d, lgtype( ncodes ), lp2d, lp3a, lp3d, lscf
      logical ls3d, lsf2d, lsf3d, ldr3d
      equivalence ( lgtype( 1 ), lp2d ), ( lgtype( 2 ), lc3d )
      equivalence ( lgtype( 3 ), lsf2d ), ( lgtype( 4 ), lp3a )
      equivalence ( lgtype( 5 ), lp3d ), ( lgtype( 6 ), ls3d )
      equivalence ( lgtype( 7 ), lsf3d ), ( lgtype( 8 ), lc2d )
      equivalence ( lgtype( 9 ), lscf ), ( lgtype( 10 ), ldr3d )
      save iuse, lgtype
c
c local variables
      integer i
      real rnew, x, y, z
      data iuse / 0 /
c
c set default for gridless codes
      ongrd1 = .true.
c change grids allowed only in hybrid mode
      if( hybrid )then
        if( iuse .eq. 0 )then
          if( ncode .gt. 9 )call crash( 'ONGRD1', 'New grid type' )
c make a local copy of the grid-type logicals and flag the smaller hybrid
          do i = 1, ncodes
            lgtype( i ) = .false.
          end do
          lgtype( igrid( 1 ) ) = .true.
          iuse = 1
        end if
c
        if( twod )then
          i = max( iz( is ), 1 )
          i = min( i, nzones )
          x = newc( 1, is ) - xcpred( 1, i, 1 )
          y = newc( 2, is ) - xcpred( 2, i, 1 )
          if( c2d )then
            ongrd1 = ( abs( x ) .lt. xm ) .and. ( abs( x ) .lt. ym )
          else
            rnew = sqrt( x * x + y * y )
            ongrd1 = rnew .lt. rgrid( 1 )
          end if
        else
          i = max( iz( is ), 1 )
          i = min( i, nzones )
          x = newc( 1, is ) - xcpred( 1, i, 1 )
          y = newc( 2, is ) - xcpred( 2, i, 1 )
          z = newc( 3, is ) - xcpred( 3, i, 1 )
          if( lc3d )then
            ongrd1 = ( abs( x ) .lt. xm ) .and.
     +               ( abs( y ) .lt. ym ) .and.
     +               ( abs( z ) .lt. zm( 1 ) )
          else if( lp3a .or. lp3d .or. lsf3d )then
            rnew = sqrt( x * x + y * y )
            ongrd1 = ( abs( z ) .lt. zm( 1 ) ) .and.
     +               ( rnew .lt. rgrid( 1 ) )
          else if( ls3d )then
            rnew = sqrt( x * x + y * y + z * z )
            ongrd1 = rnew .lt. rgrid( 1 )
          end if
        end if
      end if
      return
      end
