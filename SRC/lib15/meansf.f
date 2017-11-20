      real function meansf( rad )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the spherically averaged radial force at the given radius.  A zero
c   value is returned for a point outside the grid.  All quantities, both
c   input and returned, are in internal grid units.
c
c The force includes the contributions from rigid mass components, as well
c   as that of the active mass
c
c calling argument - radius in grid units
      real rad
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
      integer i, l, mpz
      parameter ( mpz = 30 )
      real dz, fr, fz, r, z
c
      meansf = 0
      if( ( rad .gt. 0. ) .and. ( rad .lt. rgrid( jgrid ) ) )then
c equal spacing in z - this gives equal areas on a sphere by Archimedes theorem
        dz = 2. * rad / real( mpz )
        z = .5 * dz - rad
        l = 0
        do i = 1, mpz
c skip points outside grid
          if( abs( z ) .lt. zm( jgrid ) )then
            l = l + 1
            r = sqrt( rad * rad - z * z )
c mnfrz returns the cylindrical and vertical forces averaged around a ring
            call mnfrz( r, z, fr, fz )
            meansf = meansf + ( r * fr + z * fz ) / rad
          end if
          z = z + dz
        end do
        meansf = meansf / real( l )
      end if
      return
      end

      real*8 function meansf2( rad )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the spherically averaged radial force at the given radius.  A zero
c   value is returned for a point outside the grid.  All quantities, both
c   input and returned, are in internal grid units.
c
c The force includes the contributions from rigid mass components, as well
c   as that of the active mass
c
c calling argument - radius in grid units
      real*8 rad
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
      integer i, l, mpz
      parameter ( mpz = 30 )
      real*8 dz, fr, fz, r, z
c
      meansf2 = 0
      if( ( rad .gt. 0.d0 ) .and.
     +    ( sngl( rad ) .lt. rgrid( jgrid ) ) )then
c equal spacing in z - this gives equal areas on a sphere by Archimedes theorem
        dz = 2. * rad / dble( mpz )
        z = .5 * dz - rad
        l = 0
        do i = 1, mpz
c skip points outside grid
          if( abs( z ) .lt. zm( jgrid ) )then
            l = l + 1
            r = sqrt( rad * rad - z * z )
c mnfrz returns the cylindrical and vertical forces averaged around a ring
            call mnfrz2( r, z, fr, fz )
            meansf2 = meansf2 + ( r * fr + z * fz ) / rad
          end if
          z = z + dz
        end do
        meansf2 = meansf2 / dble( l )
      end if
      return
      end
