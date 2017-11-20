      real function meanrf( rad )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Returns the azimuthally averaged radial force at the given radius in the
c   mid-plane.  A zero value is returned for a point outside the grid.  All
c   quantities, both input and returned, are in internal grid units.
c
c The force includes the contributions from rigid mass components, as well
c   as that of the active mass
c
c calling argument
      real rad
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
c external
      logical offgrd
c
c local variables
      integer m, mtp
      parameter ( mtp = 50 )
      real fz, t, x, y
      include 'inc/pi.f'
c
      if( mtp .gt. mbuff )call crash( 'MEANRF',
     +                                        'Local arrays too small' )
      newc( 1, 1 ) = rad + xcen( 1, 1 )
      newc( 2, 1 ) = xcen( 2, 1 )
      if( threed )newc( 3, 1 ) = xcen( 3, 1 )
      if( ( rad .le. 0. ) .or. offgrd( 1 ) )then
        meanrf = 0
      else
c 2-D grids
        if( twod )then
c generate a ring of test particles at this radius
          do m = 1, mtp
            t = 2. * pi * real( m - 1 ) / real( mtp )
            oldc( 1, m ) = rad * cos( t ) + xcen( 1, 1 )
            oldc( 2, m ) = rad * sin( t ) + xcen( 2, 1 )
            iz( m ) = 1
            iflag( m ) = 1
            pwt( m ) = 0
          end do
          izone = 1
          call getacc( mtp )
c average round ring
          meanrf = 0
          do m = 1, mtp
            x = oldc( 1, m ) - xcen( 1, 1 )
            y = oldc( 2, m ) - xcen( 2, 1 )
            meanrf = meanrf + ( acc( 1, m ) * x
     +                        + acc( 2, m ) * y ) / rad
          end do
c halo force already included
          meanrf = meanrf / real( mtp )
        else
c general 3-D grids
          t = 0
          call mnfrz( rad, t, meanrf, fz )
        end if
      end if
      return
      end

      real*8 function meanrf2( rad )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c Returns the azimuthally averaged radial force at the given radius in the
c   mid-plane.  A zero value is returned for a point outside the grid.  All
c   quantities, both input and returned, are in internal grid units.
c
c The force includes the contributions from rigid mass components, as well
c   as that of the active mass
c
c calling argument
      real*8 rad
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
c external
      logical offgrd
c
c local variables
      integer m, mtp
      parameter ( mtp = 50 )
      real*8 fz, t, x, y
      include 'inc/pi.f'
c
      if( mtp .gt. mbuff )call crash( 'MEANRF2',
     +                                        'Local arrays too small' )
      newc( 1, 1 ) = rad + xcen( 1, 1 )
      newc( 2, 1 ) = xcen( 2, 1 )
      if( threed )newc( 3, 1 ) = xcen( 3, 1 )
      if( ( rad .le. 0.d0 ) .or. offgrd( 1 ) )then
        meanrf2 = 0
      else
c 2-D grids
        if( twod )then
c generate a ring of test particles at this radius
          do m = 1, mtp
            t = 2. * pi * real( m - 1 ) / real( mtp )
            oldc( 1, m ) = rad * cos( t ) + xcen( 1, 1 )
            oldc( 2, m ) = rad * sin( t ) + xcen( 2, 1 )
            iz( m ) = 1
            iflag( m ) = 1
            pwt( m ) = 0
          end do
          izone = 1
          call getacc( mtp )
c average round ring
          meanrf2 = 0
          do m = 1, mtp
            x = oldc( 1, m ) - xcen( 1, 1 )
            y = oldc( 2, m ) - xcen( 2, 1 )
            meanrf2 = meanrf2 + ( acc( 1, m ) * x
     +                          + acc( 2, m ) * y ) / rad
          end do
c halo force already included
          meanrf2 = meanrf2 / real( mtp )
        else
c general 3-D grids
          t = 0
          call mnfrz2( rad, t, meanrf2, fz )
        end if
      end if
      return
      end
