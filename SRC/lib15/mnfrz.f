      subroutine mnfrz( rad, z, fr, fz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Returns the azimuthally averaged cylindrical force components at the given
c   radius and z position.  Zero values are returned for a point outside the
c   grid.  All quantities, both input and returned, are in internal grid units.
c
c With the exception of the axisymmetric polar grid, the force is evaluated
c   by creating a ring of test particles at the desired position, using
c   GETACC to determine the acceleration components for each particle and
c   averaging.  For the p3a grid, a single particle defines the average.
c
c This routine ASSUMES, for the C3D grid only, that the appropriate plane of
c   forces is ready for use
c
c calling arguments
      real fr, fz, rad, z
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
      include 'inc/model.f'
c
c externals
      logical offgrd
c
c local variables
      integer i, isave, jz, m, mtp
      parameter ( mtp = 50 )
      logical firstc, lstphv
      real t, x, y
      save firstc, isave
      include 'inc/pi.f'
      data firstc / .true. /
c
      if( firstc )then
        isave = 0
        do i = 1, ncmp
          if( igrd( i ) .eq. 1 )isave = i
        end do
        firstc = .false.
      end if
c
      if( twod )call crash( 'MNFRZ', '3-D routine called' )
c forces should be zero at centre
      if( ( rad .lt. 1.e-3 ) .and. ( abs( z ) .lt. 1.e-3 ) )then
        fr = 0
        fz = 0
      else
c set izone
        izone = 1
c no need to average for axisymmetric grid or a point on the axis
        if( p3a .or. ( rad .lt. 1.e-3 ) )then
          newc( 1, 1 ) = rad + xcen( 1, jgrid )
          newc( 2, 1 ) = xcen( 2, jgrid )
          newc( 3, 1 ) = z + xcen( 3, jgrid )
          iz( 1 ) = izone
          if( offgrd( 1 ) )iz( 1 ) = nzones
          iflag( 1 ) = isave
          call relabl( 1 )
          pwt( 1 ) = 0
c ensure point is not exactly on z-axis
          fr = max( rad, 1.e-6 )
          oldc( 1, 1 ) = fr + xcen( 1, jgrid )
          oldc( 2, 1 ) = xcen( 2, jgrid )
          oldc( 3, 1 ) = z + xcen( 3, jgrid )
          call getacc( 1 )
          fr = acc( 1, 1 )
          fz = acc( 3, 1 )
        else
c assume zone 1 so forces are not rescaled
          newc( 1, 1 ) = rad + xcen( 1, jgrid )
          newc( 2, 1 ) = xcen( 2, jgrid )
          newc( 3, 1 ) = z + xcen( 3, jgrid )
          jz = izone
          if( offgrd( 1 ) )jz = nzones
c generate a ring of test particles at the required radius and height
          do m = 1, mtp
            t = 2. * pi * real( m - 1 ) / real( mtp )
            oldc( 1, m ) = rad * cos( t ) + xcen( 1, jgrid )
            oldc( 2, m ) = rad * sin( t ) + xcen( 2, jgrid )
            oldc( 3, m ) = z + xcen( 3, jgrid )
            iz( m ) = jz
            iflag( m ) = isave
            pwt( m ) = 0
            do i = 1, 3
              newc( i, m ) = oldc( i, m )
            end do
          end do
c set hybrid label and plane number
          call relabl( mtp )
c get accelerations
          if( bht )call bhxtra( mtp )
          if( dr3d )then
            lstphv = stphev
            stphev = .false.
          end if
          call getacc( mtp )
          if( dr3d )stphev = lstphv
c average round ring
          fr = 0
          fz = 0
          do m = 1, mtp
            x = oldc( 1, m ) - xcen( 1, jgrid )
            y = oldc( 2, m ) - xcen( 2, jgrid )
            fr = fr + ( x * acc( 1, m ) + y * acc( 2, m ) )
            fz = fz + acc( 3, m )
          end do
c fixed forces already included in getacc
          fr = fr / ( rad * real( mtp ) )
          fz = fz / real( mtp )
        end if
      end if
      return
      end

      subroutine mnfrz2( rad, z, fr, fz )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c Returns the azimuthally averaged cylindrical force components at the given
c   radius and z position.  Zero values are returned for a point outside the
c   grid.  All quantities, both input and returned, are in internal grid units.
c
c With the exception of the axisymmetric polar grid, the force is evaluated
c   by creating a ring of test particles at the desired position, using
c   GETACC to determine the acceleration components for each particle and
c   averaging.  For the p3a grid, a single particle defines the average.
c
c This routine ASSUMES, for the C3D grid only, that the appropriate plane of
c   forces is ready for use
c
c calling arguments
      real*8 fr, fz, rad, z
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
      include 'inc/model.f'
c
c externals
      logical offgrd
c
c local variables
      integer i, isave, jz, m, mtp
      parameter ( mtp = 50 )
      logical firstc, lstphv
      real*8 t, x, y
      save firstc, isave
      include 'inc/pi.f'
      data firstc / .true. /
c
      if( firstc )then
        isave = 0
        do i = 1, ncmp
          if( igrd( i ) .eq. 1 )isave = i
        end do
        firstc = .false.
      end if
c
      if( twod )call crash( 'MNFRZ2', '3-D routine called' )
c forces should be zero at centre
      if( ( rad .lt. 1.d-3 ) .and. ( abs( z ) .lt. 1.d-3 ) )then
        fr = 0
        fz = 0
      else
c set izone
        izone = 1
c no need to average for axisymmetric grid or a point on the axis
        if( p3a .or. ( rad .lt. 1.d-3 ) )then
          newc( 1, 1 ) = rad + xcen( 1, jgrid )
          newc( 2, 1 ) = xcen( 2, jgrid )
          newc( 3, 1 ) = z + xcen( 3, jgrid )
          iz( 1 ) = izone
          if( offgrd( 1 ) )iz( 1 ) = nzones
          iflag( 1 ) = isave
          call relabl( 1 )
          pwt( 1 ) = 0
c ensure point is not exactly on z-axis
          fr = max( rad, 1.d-6 )
          oldc( 1, 1 ) = fr + xcen( 1, jgrid )
          oldc( 2, 1 ) = xcen( 2, jgrid )
          oldc( 3, 1 ) = z + xcen( 3, jgrid )
          call getacc( 1 )
          fr = acc( 1, 1 )
          fz = acc( 3, 1 )
        else
c assume zone 1 so forces are not rescaled
          newc( 1, 1 ) = rad + xcen( 1, jgrid )
          newc( 2, 1 ) = xcen( 2, jgrid )
          newc( 3, 1 ) = z + xcen( 3, jgrid )
          jz = izone
          if( offgrd( 1 ) )jz = nzones
c generate a ring of test particles at the required radius and height
          do m = 1, mtp
            t = 2. * pi * real( m - 1 ) / real( mtp )
            oldc( 1, m ) = rad * cos( t ) + xcen( 1, jgrid )
            oldc( 2, m ) = rad * sin( t ) + xcen( 2, jgrid )
            oldc( 3, m ) = z + xcen( 3, jgrid )
            iz( m ) = jz
            iflag( m ) = isave
            pwt( m ) = 0
            do i = 1, 3
              newc( i, m ) = oldc( i, m )
            end do
          end do
c set hybrid label and plane number
          call relabl( mtp )
c get accelerations
          if( bht )call bhxtra( mtp )
          if( dr3d )then
            lstphv = stphev
            stphev = .false.
          end if
          call getacc( mtp )
          if( dr3d )stphev = lstphv
c average round ring
          fr = 0
          fz = 0
          do m = 1, mtp
            x = oldc( 1, m ) - xcen( 1, jgrid )
            y = oldc( 2, m ) - xcen( 2, jgrid )
            fr = fr + ( x * acc( 1, m ) + y * acc( 2, m ) )
            fz = fz + acc( 3, m )
          end do
c fixed forces already included in getacc
          fr = fr / ( rad * real( mtp ) )
          fz = fz / real( mtp )
        end if
      end if
      return
      end
