      subroutine orbint( E, Lz, tol )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c subroutine to integrate an orbit of a given energy and angular momentum in
c   in a fixed, planar potential, with a specified accuracy.  The routine
c   assumes that the limiting radii and angular frequencies have already
c   been determined for the orbit.
c The three equations of motion are:
c         udot = frtot( r ) + Lz^2 / r^3
c         rdot = u
c         phidot = Lz / r^2
c   with y(1) = u, y(2) = r, y(3) = phi
c
c calling arguments
      real*8 E, Lz
      real tol
c
c common block
c
      include 'inc/orbval.f'
c
c externals
      external orbfcn, orbsto
c
c local array
      real*8 y( 3 )
c
c local variables
      integer dim, i, ifail, j, n
      logical firstc
      real*8 dt, t, tau, tol2
      parameter ( dim = 1009 )
      include 'inc/pi.f'
      save firstc
c
      data firstc / .true. /
      if( firstc )then
c store array size parameter and allocate
        if( .not. allocated( orbtab ) )then
          allocate ( orbtab( 4 * dim ) )
          idim = dim
        end if
        firstc = .false.
      end if
c
c check values of Lz and E
      if( ( Lz .ne. crrntL ) .or. ( E .ne. crrntE )
     +                 )call crash( 'ORBINT', 'Wrong values of E & Lz' )
c set time-step markers - four steps more than 1/2 period
      dt = pi / real( idim - 9 ) / omega1
      do i = 1, idim + 1
        orbtab( i ) = dt * real( i - 5 )
      end do
c total length of integration - more to ensure that last time is available
      tau = orbtab( idim ) + .5 * dt
c initial values
      is = 4
      t = 0
      y( 1 ) = 0
      y( 2 ) = rapo
      y( 3 ) = 0
c integrate motion
      n = 3
      tol2 = tol
      ifail = 0
      call ode_tab( t, tau, n, y, orbfcn, tol2, orbsto, ifail )
c set up table from before apo
      do i = 6, 9
        j = 10 - i
        orbtab( j ) = -orbtab( i )
        orbtab( idim + j ) = orbtab( idim + i )
        orbtab( 2 * idim + j ) = -orbtab( 2 * idim + i )
        orbtab( 3 * idim + j ) = -orbtab( 3 * idim + i )
      end do
      return
      end
