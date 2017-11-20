      real function meanph( rad, z )
c  Copyright (C) 2015, Jerry Sellwood
c
c The value returned by this function is the azimuthally averaged potential
c   at the given radius and z position.  The potential includes that from
c   the active mass on the grid plus contributions from any rigid mass
c   components.
c A zero value is returned for a point outside the grid.
c All quantities, both input and returned, are in internal grid units.
c
c The sign of the potential is such that its negative gradient defines the
c   gravitational acceleration.
c
c With the exception of the axisymmetric polar grid, the potential is
c   evaluated by creating a ring of test particles at the desired position,
c   using GRDPOT to determine the potential at each particle, averaging and
c   adding the fixed contributions.  For the p3a grid, a single particle
c   defines the average.
      use aarrays
      implicit none
c
c calling arguments
      real rad, z
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
c externals
      logical offgrd, ongrd1
      real axipot, grofu, halpot, phcorr, uofgr
c
c local variables
      integer j, jl, jp, jz, k, m, mtp, nring
      parameter ( mtp = 50 )
      logical lstphv
      real dr, r1, r2, t, u
      include 'inc/pi.f'
c
      meanph = 0
      if( .not. potl )call crash( 'MEANPH', 'Potentials not available' )
c set izone
      izone = 1
c 2-D polar grid
      if( p2d )then
        if( z .ne. 0. )then
          print *, 'rad, z =', rad, z
          call crash( 'MEANPH', 'Impossible calling arguments' )
        end if
        if( rad .le. rgrid( jgrid ) )then
          u = uofgr( rad )
          nring = u
          nring = max( nring + 1, 2 )
c linear interpolation in radius - not u
          r1 = grofu( real( nring - 1 ) )
          r2 = grofu( real( nring ) )
          dr = ( rad - r1 ) / ( r2 - r1 )
          j = ( nring - 1 ) * na
c average over angles
          do m = 1, na
            j = j + 1
            k = j + na
            meanph = meanph +
     +                grdfld( j, 3 ) * ( 1. - dr ) + grdfld( k, 3 ) * dr
          end do
          meanph = meanph / real( na )
c add potentials of rigid components (if any)
          if( fixrad )then
            meanph = meanph + axipot( rad )
          else
            if( suppl )meanph = meanph + phcorr( rad )
            if( rigidh )meanph = meanph + halpot( rad )
          end if
        end if
      else if( twod )then
        call crash( 'MEANPH', 'Unrecognized field method' )
c 3-D methods
      else
c 3-D axisymmetric grid
        if( p3a )then
c no need to average - field is axisymmetric
          oldc( 1, 1 ) = rad + xcen( 1, 1 )
          oldc( 2, 1 ) = xcen( 2, 1 )
          oldc( 3, 1 ) = z + xcen( 3, 1 )
          iz( 1 ) = 1
          if( offgrd( 1 ) )iz( 1 ) = nzones
          label( 1 ) = 1
          if( hybrid )then
            newc( 1, 1 ) = oldc( 1, 1 )
            newc( 2, 1 ) = oldc( 2, 1 )
            newc( 3, 1 ) = oldc( 3, 1 )
            if( .not. ongrd1( 1 ) )label( 1 ) = 2
          end if
          iflag( 1 ) = 1
          pwt( 1 ) = 0
          call getacc( 1 )
          meanph = gpot( 1 )
        else
c fully 3-D grids
          if( mtp .gt. mbuff )call crash( 'MEANPH',
     +                                       'Local arrays too small' )
          newc( 1, 1 ) = rad + xcen( 1, 1 )
          newc( 2, 1 ) = xcen( 2, 1 )
          newc( 3, 1 ) = z + xcen( 3, 1 )
          jz = 1
          if( offgrd( 1 ) )jz = nzones
          jl = 1
          if( hybrid .and. ( .not. ongrd1( 1 ) ) )jl = 2
          jp = 1
          if( nplanes .gt. 1 )jp = z + zm( jgrid ) + 1.
c no need to average for a point on the axis
          k = mpt
          if( rad .eq. 0. )k = 1
c generate a ring of test particles at specified z height
          do m = 1, k
            t = 2. * pi * real( m - 1 ) / real( k )
            oldc( 1, m ) = rad * cos( t ) + xcen( 1, jl )
            oldc( 2, m ) = rad * sin( t ) + xcen( 2, jl )
            oldc( 3, m ) = z + xcen( 3, jl )
            iz( m ) = jz
            iflag( m ) = 1
            label( m ) = jl
            pwt( m ) = 0
            ipln( m ) = jp
          end do
          if( bht )call bhxtra( k )
          if( dr3d )then
            lstphv = stphev
            stphev = .false.
          end if
          call getacc( k )
          if( dr3d )stphev = lstphv
c average round ring
          do m = 1, k
            meanph = meanph + gpot( m )
          end do
          meanph = meanph / real( k )
        end if
c add potentials of rigid components (if any)
        if( fixrad )then
          meanph = meanph + axipot( rad )
        else
          r1 = sqrt( rad * rad + z * z )
          if( suppl )meanph = meanph + phcorr( r1 )
          if( rigidh )meanph = meanph + halpot( r1 )
        end if
      end if
      return
      end

      real*8 function meanph2( rad, z )
c  Copyright (C) 2015, Jerry Sellwood
c
c The value returned by this function is the azimuthally averaged potential
c   at the given radius and z position.  The potential includes that from
c   the active mass on the grid plus contributions from any rigid mass
c   components.
c A zero value is returned for a point outside the grid.
c All quantities, both input and returned, are in internal grid units.
c
c The sign of the potential is such that its negative gradient defines the
c   gravitational acceleration.
c
c With the exception of the axisymmetric polar grid, the potential is
c   evaluated by creating a ring of test particles at the desired position,
c   using GRDPOT to determine the potential at each particle, averaging and
c   adding the fixed contributions.  For the p3a grid, a single particle
c   defines the average.
      use aarrays
      implicit none
c
c calling arguments
      real*8 rad, z
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
c externals
      logical offgrd, ongrd1
      real axipot, grofu, halpot, phcorr, uofgr
c
c local variables
      integer j, jl, jp, jz, k, m, mtp, nring
      parameter ( mtp = 50 )
      logical lstphv
      real t, u
      real*8 dr, r1, r2
      include 'inc/pi.f'
c
      meanph2 = 0
      if( .not.potl )call crash( 'MEANPH2', 'Potentials not available' )
c set izone
      izone = 1
c 2-D polar grid
      if( p2d )then
        if( z .ne. 0.d0 )then
          print *, 'rad, z =', rad, z
          call crash( 'MEANPH2', 'Impossible calling arguments' )
        end if
        if( sngl( rad ) .le. rgrid( jgrid ) )then
          u = uofgr( sngl( rad ) )
          nring = u
          nring = max( nring + 1, 2 )
c linear interpolation in radius - not u
          r1 = grofu( real( nring - 1 ) )
          r2 = grofu( real( nring ) )
          dr = ( rad - r1 ) / ( r2 - r1 )
          j = ( nring - 1 ) * na
c average over angles
          do m = 1, na
            j = j + 1
            k = j + na
            meanph2 = meanph2 +
     +              grdfld( j, 3 ) * ( 1.d0 - dr ) + grdfld( k, 3 ) * dr
          end do
          meanph2 = meanph2 / real( na )
c add potentials of rigid components (if any)
          if( fixrad )then
            meanph2 = meanph2 + axipot( sngl( rad ) )
          else
            if( suppl )meanph2 = meanph2 + phcorr( sngl( rad ) )
            if( rigidh )meanph2 = meanph2 + halpot( sngl( rad ) )
          end if
        end if
      else if( twod )then
        call crash( 'MEANPH2', 'Unrecognized field method' )
c 3-D methods
      else
c 3-D axisymmetric grid
        if( p3a )then
c no need to average - field is axisymmetric
          oldc( 1, 1 ) = rad + xcen( 1, 1 )
          oldc( 2, 1 ) = xcen( 2, 1 )
          oldc( 3, 1 ) = z + xcen( 3, 1 )
          iz( 1 ) = 1
          if( offgrd( 1 ) )iz( 1 ) = nzones
          label( 1 ) = 1
          if( hybrid )then
            newc( 1, 1 ) = oldc( 1, 1 )
            newc( 2, 1 ) = oldc( 2, 1 )
            newc( 3, 1 ) = oldc( 3, 1 )
            if( .not. ongrd1( 1 ) )label( 1 ) = 2
          end if
          iflag( 1 ) = 1
          pwt( 1 ) = 0
          call getacc( 1 )
          meanph2 = gpot( 1 )
        else
c fully 3-D grids
          if( mtp .gt. mbuff )call crash( 'MEANPH2',
     +                                       'Local arrays too small' )
          newc( 1, 1 ) = rad + xcen( 1, 1 )
          newc( 2, 1 ) = xcen( 2, 1 )
          newc( 3, 1 ) = z + xcen( 3, 1 )
          jz = 1
          if( offgrd( 1 ) )jz = nzones
          jl = 1
          if( hybrid .and. ( .not. ongrd1( 1 ) ) )jl = 2
          jp = 1
          if( nplanes .gt. 1 )jp = z + zm( jgrid ) + 1.
c no need to average for a point on the axis
          k = mpt
          if( rad .eq. 0.d0 )k = 1
c generate a ring of test particles at specified z height
          do m = 1, k
            t = 2. * pi * real( m - 1 ) / real( k )
            oldc( 1, m ) = rad * cos( t ) + xcen( 1, jl )
            oldc( 2, m ) = rad * sin( t ) + xcen( 2, jl )
            oldc( 3, m ) = z + xcen( 3, jl )
            iz( m ) = jz
            iflag( m ) = 1
            label( m ) = jl
            pwt( m ) = 0
            ipln( m ) = jp
          end do
          if( bht )call bhxtra( k )
          if( dr3d )then
            lstphv = stphev
            stphev = .false.
          end if
          call getacc( k )
          if( dr3d )stphev = lstphv
c average round ring
          do m = 1, k
            meanph2 = meanph2 + gpot( m )
          end do
          meanph2 = meanph2 / real( k )
        end if
c add potentials of rigid components (if any)
        if( fixrad )then
          meanph2 = meanph2 + axipot( sngl( rad ) )
        else
          r1 = sqrt( rad * rad + z * z )
          if( suppl )meanph2 = meanph2 + phcorr( sngl( r1 ) )
          if( rigidh )meanph2 = meanph2 + halpot( sngl( r1 ) )
        end if
      end if
      return
      end
