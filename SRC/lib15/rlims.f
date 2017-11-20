      subroutine rlims( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns values of apo and peri in z=0 plane for a given E and Lz
c
c calling arguments
      real*8 E, Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      include 'inc/orbval.f'
c
c externals
      real*8 akappa, phitot, rcirc
c
c local arrays
      real*8 w( 4 )
c
c local variables
      integer ierr, ifail, ind, ir, irad, jrad
      real*8 ak, Ec, r, rc, r2, tol, tiny, u2
      parameter ( tiny = 1.d-10 )
      save ierr
      data ierr / 0 /
c
c store new values
      crrntE = E
      crrntL = Lz
c start from circular orbit
      rc = rcirc( Lz )
      Ec = phitot( rc )
      if( rc .gt. 0.d0 )Ec = Ec + .5 * ( Lz / rc )**2
      ak = akappa( rc )
c check that a solution exists
      if( E .le. Ec )then
        rperi = rc
        rapo = rc
        if( ( rc .gt. tiny ) .and. abs( E - Ec ) .gt. .000001d0 )then
          ierr = ierr + 1
          print *, 'in rlims E, Lz =', E, Lz
          print *, ierr, rc, abs( E - Ec )
          if( ierr .gt. 20 )call crash( 'RLIMS', 'Impossible E and Lz' )
        end if
      else if( E - Ec .lt. .5d-8 * ( rc * ak )**2 )then
c epicycle approx
        r2 = sqrt( 2. * ( E - Ec ) ) / ak
        rperi = rc - r2
        rapo = rc + r2
      else
c zero angular momentum limit
        if( Lz .eq. 0. )then
          rperi = 0.
          jrad = 2
        else
          jrad = 1
        end if
c find both extrema
        do irad = jrad, 2
c choose radii to bracket the zero
          r = rc
          r2 = rc
          if( rc .eq. 0. )then
            r = 1.e-20
            r2 = 1
          end if
    1     if( irad .eq. 1 )r2 = .09999 * r2
          if( irad .eq. 2 )r2 = 10.001 * r2
          if( ( ncmp .eq. 1 ) .and. ( ( ctype( 1 ) .eq. 'MFK ' )
     +           .or. ( ctype( 1 ) .eq. 'SPLY' ) ) )r2 = min( r2, rmax )
c ask for impossibly high precision
          tol = 1.d-30
          ir = 1
          ind = 1
          ifail = 1
c find zero
          do while ( ind .gt. 0 )
            call fndzro( r, r2, u2, tol, ir, w, ind, ifail )
            u2 = 2. * ( E - phitot( r ) ) - ( Lz / r )**2
          end do
c check ifail
          if( ifail .eq. 1 )then
            u2 = max( r, r2 )
            r = min( r, r2 )
            if( irad .eq. 1 )then
              r2 = r
              r = ( 1. + tiny ) * u2
            else
              r = ( 1. + tiny ) * r
              r2 = u2
            end if
c try wider radial range unless zero is very close to rmax
            if( ( abs( r2 - rmax ) .gt. 1.d-6 ) .and.
     +          ( r2 .lt. 100.d0 * rmax ) )go to 1
            u2 = 2. * ( E - phitot( r2 ) ) - ( Lz / r2 )**2
            if( abs( u2 ) .gt. 1.d-5 )then
              print *, 'Lz, E =', Lz, E
              print *, r, 2. * ( E - phitot( r ) ) - ( Lz / r )**2
              print *, r2, 2. * ( E - phitot( r2 ) ) - ( Lz / r2 )**2
              call crash( 'RLIMS', 'No solution found' )
            end if
            r = rmax
          end if
c store result
          if( irad .eq. 1 )rperi = r
          if( irad .eq. 2 )rapo = r
        end do
      end if
      return
      end
