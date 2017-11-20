      logical function offgrd( is )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Determines whether particle is outside the grid on the basis of its newc
c   position
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
c local "common"
      integer iuse, jg
      logical lc2d, lc3d, lgtype( 0:ncodes ), lp2d, lp3a, lp3d, lscf
      logical ls3d, lsf2d, lsf3d, ldr3d, lbht, lnone
      equivalence ( lgtype( 1 ), lp2d ), ( lgtype( 2 ), lc3d )
      equivalence ( lgtype( 3 ), lsf2d ), ( lgtype( 4 ), lp3a )
      equivalence ( lgtype( 5 ), lp3d ), ( lgtype( 6 ), ls3d )
      equivalence ( lgtype( 7 ), lsf3d ), ( lgtype( 8 ), lc2d )
      equivalence ( lgtype( 9 ), lscf ), ( lgtype( 10 ), ldr3d )
      equivalence ( lgtype( 11 ), lbht ), ( lgtype( 0 ), lnone )
      save iuse, jg, lgtype
c
c local variables
      integer i
      real rnew, x, y, z
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
        if( ncode .gt. 11 )call crash( 'OFFGRD',
     +                                        'Unrecognised grid type' )
c make a local copy of the grid-type logicals and select the largest
        do i = 0, ncodes
          lgtype( i ) = .false.
        end do
c largest grid is grid 1 for twogrd and ngrid for hybrid
        jg = 1
        if( hybrid )jg = ngrid
        lgtype( igrid( jg ) ) = .true.
        iuse = 1
      end if
c set default for gridless codes
      offgrd = .false.
      if( twod )then
        x = newc( 1, is ) - xcpred( 1, nzones, jg )
        y = newc( 2, is ) - xcpred( 2, nzones, jg )
        if( c2d )then
          offgrd = ( abs( x ) .gt. xm ) .or. ( abs( y ) .gt. ym )
        else
          rnew = sqrt( x * x + y * y )
          offgrd = rnew .gt. rgrid( jg )
        end if
      else
        x = newc( 1, is ) - xcpred( 1, nzones, jg )
        y = newc( 2, is ) - xcpred( 2, nzones, jg )
        z = newc( 3, is ) - xcpred( 3, nzones, jg )
        if( lc3d )then
          offgrd = ( abs( x ) .gt. xm ) .or.
     +             ( abs( y ) .gt. ym ) .or.
     +             ( abs( z ) .gt. zm( jg ) )
        else if( lp3a .or. lp3d .or. lsf3d )then
          rnew = sqrt( x * x + y * y )
          offgrd = ( abs( z ) .gt. zm( jg ) )
     +             .or. ( rnew .gt. rgrid( jg ) )
        else if( ls3d )then
          rnew = sqrt( x * x + y * y + z * z )
          offgrd = rnew .gt. rgrid( jg )
        end if
      end if
      return
      end
