      subroutine rhoprf
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to calculate and save the radial density profile of an assumed
c  spherical model.  It sorts the particles in radius and calculates the
c  mean local density from the radial distance between fixed numbers of
c  particles
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c local allocatable arrays
      integer, allocatable :: iprf(:)
c
      real, allocatable :: prfrad(:)
c
      real, allocatable :: w(:)
c
c local variables
      character*4 bstr( 2 )
      integer i, ia, ip, j, jp, lw, n
      real a, ahol, fac, r1, r2, x, y, z
      equivalence ( a, ia )
      include 'inc/pi.f'
c
      data bstr / 'SIGR', 'RHOR' /
c
      if( parallel )call crash( 'RHOPRF', 'mpi version needed' )
c allocate space
      n = 0
      do ip = 1, ncmp
        n = max( n, nsp( ip ) )
      end do
      allocate ( iprf( n ) )
      allocate ( prfrad( n ) )
      lw = 2 * ( ( n - 1 ) / mprho + 1 )
c work over active populations
      do ip = 1, ncmp
        n = 0
        do ilist = 1, nlists
          i = islist( 1, ilist, myid + 1 )
          do while ( i .ge. 0 )
c get population flag
            a = ptcls( i + nwpp - 1 )
            if( ia .eq. ip )then
              n = n + 1
c get squared radii of all particles
              if( disc( ip ) )then
                x = ptcls( i + 1 ) - xcen( 1, ip )
                y = ptcls( i + 2 ) - xcen( 2, ip )
c may need to time center positions
                if( .not. lfrv )then
                  j = ncoor / 2
                  x = x - .5 * ptcls( i + j + 1 )
                  y = y - .5 * ptcls( i + j + 2 )
                end if
                prfrad( n ) = x**2 + y**2
              else
                x = ptcls( i + 1 ) - xcen( 1, ip )
                y = ptcls( i + 2 ) - xcen( 2, ip )
                z = ptcls( i + 3 ) - xcen( 3, ip )
c may need to time center positions
                if( .not. lfrv )then
                  x = x - .5 * ptcls( i + 4 )
                  y = y - .5 * ptcls( i + 5 )
                  z = z - .5 * ptcls( i + 6 )
                end if
                prfrad( n ) = x**2 + y**2 + z**2
              end if
            end if
            a = ptcls( i + nwpp )
            i = ia
          end do
        end do
c
        if( n .gt. mprho )then
          allocate ( w( lw ) )
          if( master )then
c rank particles in radius
            call rnkmrg( prfrad, n, iprf )
            j = 2 * ( n / mprho )
c bin up the results
            do i = 1, j
              w( i ) = 0
            end do
            j = 0
            do i = 1, n
              w( j + 2 ) = w( j + 2 ) + 1
c save radius of selected particles
              if( mod( i, mprho ) .eq. 0 )then
                jp = iprf( i )
                w( j + 1 ) = sqrt( prfrad( jp ) ) / lscale
                j = j + 2
              end if
            end do
          end if
          j = 2 * ( n / mprho )
c may need to take individual particle masses into account
          if( uqmass )then
c copy masses of particles into the work array in the same order as radii
            n = 0
            do ilist = 1, nlists
              i = islist( 1, ilist, myid + 1 )
              do while ( i .ge. 0 )
                a = ptcls( i + nwpp - 1 )
                if( ia .eq. ip )then
                  n = n + 1
                  prfrad( n ) = ptcls( i + ncoor + 1 )
                end if
                a = ptcls( i + nwpp )
                i = ia
              end do
            end do
c revise masses in each shell
            do i = 2, j, 2
              w( i ) = 0
            end do
            if( master)then
              j = 0
              do i = 1, n
                jp = iprf( i )
                w( j + 2 ) = w( j + 2 ) + prfrad( jp )
                if( mod( i, mprho ) .eq. 0 )j = j + 2
              end do
            end if
          end if
c fac is the mass scale for one particle in external units
          fac = pmass / ( lscale**3 * ts**2 )
          if( master )then
c compute density profiles
            j = 2 * ( n / mprho )
            r2 = 0
            if( disc( ip ) )then
c mass per unit area
              fac = fac / pi
              do i = 1, j, 2
                r1 = r2
                r2 = w( i )
                w( i ) = .5 * ( r1 + r2 )
                w( i + 1 ) = fac * w( i + 1 ) / ( r2**2 - r1**2 )
              end do
              read( bstr( 1 ), '( a4 )' )ahol
              write( nphys )ip, ahol, istep, j, mprho
              write( nphys )( w( i ), i = 1, j )
            else
c mass per unit volume
              fac = 3. * fac / ( 4. * pi )
              do i = 1, j, 2
                r1 = r2
                r2 = w( i )
                w( i ) = .5 * ( r1 + r2 )
                w( i + 1 ) = fac * w( i + 1 ) / ( r2**3 - r1**3 )
              end do
              read( bstr( 2 ), '( a4 )' )ahol
              write( nphys )ip, ahol, istep, j, mprho
              write( nphys )( w( i ), i = 1, j )
            end if
          end if
          deallocate ( w )
        end if
      end do
      deallocate ( iprf )
      deallocate ( prfrad )
      return
      end
