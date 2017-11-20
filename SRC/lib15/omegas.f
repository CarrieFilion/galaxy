      subroutine omegas( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the rate of change of the angle variables for a given E & Lz
c   uses Burkardt's routine JACOBI_EK_COMPUTE to replace NAG D01BCF
c
c calling arguments
      real*8 E, Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/orbval.f'
c
c externals
      real*8 akappa, Lzmax, omegac, omgsfn, quad_gnr, rcirc, taufn
      external taufn
c
c local arrays
      integer npts
      parameter ( npts = 128 )
      real*8 abscis( npts ), weight( npts )
      save abscis, weight
c
c local variables
      integer i, ier, ifail
      logical firstc
      real*8 angle, eccen, epsa, epsr, nghalf, Lz1, r, tau
      include 'inc/pi.f'
      parameter ( nghalf = -0.5d0 )
      save firstc
c
      data firstc / .true. /
c
c find peri and apo
      call rlims( E, Lz )
      if( rapo + rperi .le. 0.d0 )then
        eccen = 0
      else
        eccen = ( rapo - rperi ) / ( rapo + rperi )
      end if
c use epicycle approximation for nearly circular orbits
      if( eccen .lt. .001 )then
c nearly circular orbits
        Lz1 = Lzmax( E )
        r = rcirc( Lz1 )
        omega1 = akappa( r )
        omega2 = omegac( r )
        omega2 = sign( omega2, Lz )
      else
        if( crrntL .eq. 0.d0 )then
c special treatment for radial orbits
          epsa = 1.d-8
          epsr = epsa
          ier = 1
          tau = quad_gnr( taufn, rperi, rapo, epsa, epsr, ier )
          if( ier .gt. 0 )then
            print *, 'ier =', ier, ' from QUAD_GNR'
            call crash( 'OMEGAS', 'QUADPACK error' )
          end if
c radial orbits
          angle = 0
        else
c set weights and abscissae for non-adaptive Gauss-Jacobi quadrature
c   for (x-a)**alfa*(b-x)**beta
          if( firstc )then
c quadrature rule for integration from -1 to 1
            call
     +         jacobi_ek_compute( npts, nghalf, nghalf, abscis, weight )
            firstc = .false.
          end if
c evaluate sum
          tau = 0
          angle = 0
          do i = 1, npts
c rescale abscissae only for integration from a to b as
c     x(a,b) = 0.5 * (b-a) * x(-1,+1) + 0.5 * (a+b)
            r = .5d0 * ( rapo - rperi ) * abscis( i ) +
     +          .5d0 * ( rapo + rperi )
            tau = tau + weight( i ) * taufn( r )
            angle = angle + weight( i ) * omgsfn( r )
          end do
        end if
        if( Lz .ge. 0. )angle = max( angle, .5 * pi )
        if( Lz .lt. 0. )angle = min( angle, -.5 * pi )
c tau is half radial period
        omega1 = pi / tau
c both angle and tau are computed over half a radial oscillation
        omega2 = angle / tau
      end if
      return
      end
