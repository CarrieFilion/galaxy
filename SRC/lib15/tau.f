      real*8 function tau( E, Lz )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c returns the radial period for a particle in the z=0 plane
c    for a given E & Lz
c uses Burkardt's routine JACOBI_EK_COMPUTE to replace NAG D01BCF
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
      real*8 akappa, Lzmax, quad_gnr, rcirc, taufn
      external taufn
c
c local arrays
      integer npts
      parameter ( npts = 32 )
      real*8 abscis( npts ), weight( npts )
      save abscis, weight
c
c local variables
      integer i, ier
      logical firstc
      real*8 eccen, epsa, epsr, Lz1, nghalf, r
      save firstc
      parameter ( nghalf = -0.5d0 )
      include 'inc/pi.f'
c
      data firstc / .true. /
c
c find apo and peri
      call rlims( E, Lz )
      crrntE = E
      crrntL = Lz
c special case of particle at rest at centre
      eccen = 0
      if( rapo .gt. 0. )eccen = ( rapo - rperi ) / ( rapo + rperi )
      if( eccen .gt. .001 )then
c accurate radial period
        if( crrntL .eq. 0.d0 )then
c special treatment for radial orbits
          epsa = 1.d-8
          epsr = epsa
          ier = 1
          tau = quad_gnr( taufn, rperi, rapo, epsa, epsr, ier )
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
          do i = 1, npts
c rescale abscissae only for integration from a to b as
c     x(a,b) = 0.5 * (b-a) * x(-1,+1) + 0.5 * (a+b)
            r = .5d0 * ( rapo - rperi ) * abscis( i ) +
     +          .5d0 * ( rapo + rperi )
            tau = tau + weight( i ) * taufn( r )
          end do
        end if
        tau = 2. * tau
      else
c epicyclic period will do
        Lz1 = Lzmax( E )
        r = rcirc( Lz1 )
        tau = 2. * pi / akappa( r )
      end if
      return
      end
